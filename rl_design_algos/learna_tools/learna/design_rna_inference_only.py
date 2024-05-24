import time
import warnings
warnings.filterwarnings('ignore')

import os
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'  # '2'

from tqdm import tqdm

import tensorflow as tf
from learna_tools.tensorforce.runner import Runner

from learna_tools.learna.agent import NetworkConfig, get_network, AgentConfig, get_agent_fn
from learna_tools.learna.environment import RnaDesignEnvironment, RnaDesignEnvironmentConfig

from concurrent.futures import ThreadPoolExecutor, as_completed

SOLUTIONS = []
CANDIDATES = {}


def _get_episode_finished(timeout, stop_once_solved, num_solutions, hamming_tolerance, pbar=None, solutions=None):  # , SOLUTIONS):
    """
    Check for timeout after each episode of designing one entire target structure.

    Args:
        timeout: Maximum time allowed to solve one target structure.
        stop_once_solved: Defines if agent should stop after solving a target structure.

    Returns:
        episode_finish: Inner function that handles timeout.
    """
    start_time = time.time()

    def episode_finished(runner):
        env = runner.environment
        target_id = env.episodes_info[-1].target_id
        candidate_solution = env.design.primary
        last_reward = runner.episode_rewards[-1]
        last_fractional_hamming = env.episodes_info[-1].normalized_hamming_distance
        elapsed_time = time.time() - start_time
        hamming_distance = env.episodes_info[-1].hamming_distance
        structure = env.episodes_info[-1].structure
        if solutions is None:
            global CANDIDATES
            global SOLUTIONS
        else:
            SOLUTIONS = solutions

        if last_reward == 1.0 or hamming_distance <= hamming_tolerance:
            # use only unique solutions
            if candidate_solution not in [s['sequence'] for s in SOLUTIONS]:
                # CANDIDATES[candidate_solution] = True
                SOLUTIONS.append({'Id': target_id,
                                  'time': elapsed_time,
                                  'hamming_distance': hamming_distance,
                                  'rel_hamming_distance': last_fractional_hamming,
                                  'sequence': candidate_solution,
                                  'structure': structure,
                                  })
                if pbar is not None:
                    pbar.update(1)

        # print(elapsed_time, last_reward, last_fractional_hamming, candidate_solution)
        num_solutions_satisfied = len(SOLUTIONS) >= num_solutions
        # print(num_solutions_satisfied)
        no_timeout = not timeout or elapsed_time < timeout
        stop_since_solved = stop_once_solved and last_reward == 1.0
        keep_running = not stop_since_solved and no_timeout and not num_solutions_satisfied
        # print(keep_running, stop_since_solved, no_timeout, num_solutions_satisfied)
        return keep_running
    return episode_finished

def worker(config, target):
    solutions = []
    environment = RnaDesignEnvironment([target], config.env_config)
    get_agent = get_agent_fn(
        environment=environment,
        network=get_network(config.network_config),
        agent_config=config.agent_config,
        
        session_config=tf.ConfigProto(
            intra_op_parallelism_threads=1,
            inter_op_parallelism_threads=1,
            allow_soft_placement=True,
            # device_count={"CPU": 1},
        ),
        
        restore_path=config.restore_path,
    )
    runner = Runner(get_agent, environment)
    runner.run(
        deterministic=False,
        restart_timeout=config.restart_timeout,
        stop_learning=config.stop_learning,
        episode_finished=_get_episode_finished(config.timeout,
                                               False,
                                               config.num_solutions,
                                               config.hamming_tolerance,
                                               solutions,
                                               # SOLUTIONS,
                                               ),
    )
    return SOLUTIONS

def design_rna(
    dot_brackets,
    timeout,
    restore_path,
    stop_learning,
    restart_timeout,
    network_config,
    agent_config,
    env_config,
    num_solutions,
    hamming_tolerance=0,
    share_agent=True,
    num_cores=1,
):
    """
    Main function for RNA design. Instantiate an environment and an agent to run in a
    tensorforce runner.

    Args:
        TODO
        timeout: Maximum time to run.
        restore_path: Path to restore saved configurations/models from.
        stop_learning: If set, no weight updates are performed (Meta-LEARNA).
        restart_timeout: Time interval for restarting of the agent.
        network_config: The configuration of the network.
        agent_config: The configuration of the agent.
        env_config: The configuration of the environment.

    Returns:
        Episode information.
    """
    print('### Start Predictions')
    print()
    print('### Solutions per target:', num_solutions)
    print('### Number of targets:', len(dot_brackets))
    print('### With Ids:', ', '.join([str(i) for i, _ in dot_brackets]))
    print('### Timeout per target:', timeout)
    print()
    
    if num_cores > 1:
        config = {
            'dot_brackets': dot_brackets,
            'timeout': timeout,
            'restore_path': restore_path,
            'stop_learning': stop_learning,
            'restart_timeout': restart_timeout,
            'network_config': network_config,
            'agent_config': agent_config,
            'env_config': env_config,
            'num_solutions': num_solutions,
            'hamming_tolerance': hamming_tolerance,
            'share_agent': share_agent
        }
        
        with ThreadPoolExecutor(max_workers=num_cores) as executor:
            futures = {executor.submit(worker, config, target=target): target for target in dot_brackets}
            for future in as_completed(futures):
                result = future.result()
                yield futures[future], result

    else:
    
        session_config = tf.ConfigProto(
                intra_op_parallelism_threads=1,
                inter_op_parallelism_threads=1,
                allow_soft_placement=True,
                device_count={"CPU": 1},
            )
        
        env_config.use_conv = any(map(lambda x: x > 1, network_config.conv_sizes))
        env_config.use_embedding = bool(network_config.embedding_size)
        network = get_network(network_config)
        
        if share_agent:
            dummy_env = RnaDesignEnvironment([dot_brackets[0]], env_config)    
        
            get_agent = get_agent_fn(
                environment=dummy_env,
                network=network,
                agent_config=agent_config,
                session_config=session_config,
                restore_path=restore_path,
            )
    
        
    
        for i, target in dot_brackets:
            print(f'### Process target {i}')
            print()
            
            pbar = tqdm(total=num_solutions)
        
            environment = RnaDesignEnvironment([(i, target)], env_config)
    
            if not share_agent:
            
                get_agent = get_agent_fn(
                    environment=environment,
                    network=network,
                    agent_config=agent_config,
                    session_config=session_config,
                    restore_path=restore_path,
                )
        
            global SOLUTIONS
        
            SOLUTIONS = []
        
            # Runner restarts the agent by calling get_agent again
            
            
            runner = Runner(get_agent, environment)
        
            # stop_once_solved = len(dot_brackets) == 1
            stop_once_solved = False
            runner.run(
                deterministic=False,
                restart_timeout=restart_timeout,
                stop_learning=stop_learning,
                episode_finished=_get_episode_finished(timeout,
                                                       stop_once_solved,
                                                       num_solutions,
                                                       hamming_tolerance,
                                                       pbar,
                                                       # SOLUTIONS,
                                                       ),
            )
            pbar.close()
            if len(SOLUTIONS) < num_solutions:
                print(f'WARNING: Found {len(SOLUTIONS)} solutions for target {i}' +
                      f' but wanted {num_solutions}.')
                print('Please consider increasing the timeout or decreasing the number of solutions' + "\033[0m.")
            yield i, SOLUTIONS
        
    
if __name__ == "__main__":
    import argparse
    from pathlib import Path
    from learna_tools.learna.data.parse_dot_brackets import parse_dot_brackets

    parser = argparse.ArgumentParser()

    # Data
    parser.add_argument(
        "--target_structure_path", type=Path, help="Path to sequence to run on"
    )
    parser.add_argument("--data_dir", default="data", help="Data directory")
    parser.add_argument("--dataset", type=Path, help="Available: eterna, rfam_taneda")
    parser.add_argument(
        "--target_structure_ids",
        default=None,
        required=False,
        type=int,
        nargs="+",
        help="List of target structure ids to run on",
    )

    # Model
    parser.add_argument("--restore_path", type=Path, help="From where to load model")
    parser.add_argument("--stop_learning", action="store_true", help="Stop learning")
    parser.add_argument("--random_agent", action="store_true", help="Use random agent")

    # Timeout behaviour
    parser.add_argument("--timeout", default=None, type=int, help="Maximum time to run")

    # Hyperparameters
    parser.add_argument("--learning_rate", type=float, help="Learning rate to use")
    parser.add_argument(
        "--mutation_threshold", type=int, help="Enable MUTATION with set threshold"
    )
    parser.add_argument(
        "--reward_exponent", default=1, type=float, help="Exponent for reward shaping"
    )
    parser.add_argument(
        "--state_radius", default=0, type=int, help="Radius around current site"
    )
    parser.add_argument(
        "--conv_sizes", type=int, default=[1], nargs="+", help="Size of conv kernels"
    )
    parser.add_argument(
        "--conv_channels", type=int, default=[50], nargs="+", help="Channel size of conv"
    )
    parser.add_argument(
        "--num_fc_layers", type=int, default=2, help="Number of FC layers to use"
    )
    parser.add_argument(
        "--fc_units", type=int, default=50, help="Number of units to use per FC layer"
    )
    parser.add_argument(
        "--batch_size", type=int, default=8, help="Batch size for ppo agent"
    )
    parser.add_argument(
        "--entropy_regularization", type=float, default=1.5e-3, help="The output entropy"
    )
    parser.add_argument(
        "--restart_timeout", type=int, help="Time after which to restart the agent"
    )
    parser.add_argument("--lstm_units", type=int, help="The number of lstm units")
    parser.add_argument("--num_lstm_layers", type=int, help="The number of lstm layers")
    parser.add_argument("--embedding_size", type=int, help="The size of the embedding")

    parser.add_argument("--min_solutions", type=int, default=1, help="Number of optimal solutions")

    args = parser.parse_args()

    network_config = NetworkConfig(
        conv_sizes=args.conv_sizes,  # radius * 2 + 1
        conv_channels=args.conv_channels,
        num_fc_layers=args.num_fc_layers,
        fc_units=args.fc_units,
        lstm_units=args.lstm_units,
        num_lstm_layers=args.num_lstm_layers,
        embedding_size=args.embedding_size,
    )
    agent_config = AgentConfig(
        learning_rate=args.learning_rate,
        batch_size=args.batch_size,
        entropy_regularization=args.entropy_regularization,
        random_agent=args.random_agent,
    )
    env_config = RnaDesignEnvironmentConfig(
        mutation_threshold=args.mutation_threshold,
        reward_exponent=args.reward_exponent,
        state_radius=args.state_radius,
    )
    dot_brackets = parse_dot_brackets(
        dataset=args.dataset,
        data_dir=args.data_dir,
        target_structure_ids=args.target_structure_ids,
        target_structure_path=args.target_structure_path,
    )

    design_rna(
        dot_brackets,
        timeout=args.timeout,
        restore_path=args.restore_path,
        stop_learning=args.stop_learning,
        restart_timeout=args.restart_timeout,
        network_config=network_config,
        agent_config=agent_config,
        env_config=env_config,
        num_solutions=args.min_solutions,
        hamming_tolerance=0,
        share_agent=True,
    )
