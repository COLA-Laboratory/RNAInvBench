<h1>RNAInvBench: Benchmarking RNA Inverse Folding with Secondary and
Tertiary Structures</h1>
<p>RNAInvBench is a Python framework that provides an environment, data, and baseline algorithms for training/testing of the RNA Inverse Design Problem.
The package contains the following algorithms:

- antaRNA
- aRNAque
- IncaRNAtion
- LEARNA
- Meta-LEARNA
- Meta-LEARNA-Adapt
- RNAInverse
- SAMFEO
- SentRNA

This package has been set-out for solving for sequences that are pseudoknot-free and contain pseudoknots. We use the pKiss package for our folding engine on all pseudoknot related tasks, and we use ViennaRNA 2 for our folding engine on all pseudoknot-free related tasks.

The package contains the following datasets to be used for training and testing:
- Eterna100-v1
- Eterna100-v2
- Pseudobase++
- Rfam-Learn-Test
- Rfam-Learn-Train
- Rfam-Taneda
- Rfam-Test
- RNA-strand

The package contains dockerfiles for each environment that is required to use the algorithms. Please make use of these when using the package!

We provide example code to run the algorithms, each example run uses the following naming convention: [algorithm_name]_run.py

</p>

<h2>Installation Instructions</h2>

<p>1. Install RNAInvBench source files.</p>
<p>2. Install <a href="https://docs.docker.com/get-docker/">docker</a> (if not done so already).</p>
<p>3. Navigate to the project directory.</p>
<p>4. Run <b>docker-compose up --build -d</b> to build the docker container.</p>
Once the docker-container has been fully built, you should see a new container in your docker application. You can then use the following commands to run the example main.py scripts provided:

<p>NOTE: RL represents the particular docker container to use, you will need to choose between rl, optim and antarna.</p>
<p>NOTE2: rnainvbench-main-rl-1 represents the name of the RL docker container, your name may differ - you should check the docker application for your own unique name.</p>
<p><b>docker exec -it $(docker-compose ps -q rl) rnainvbench-main-rl-1 python /app/rl_design_algos/main.py</b></p>
