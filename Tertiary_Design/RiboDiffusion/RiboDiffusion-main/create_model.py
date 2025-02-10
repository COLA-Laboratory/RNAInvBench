from models.utils import create_model
"""Training"""
import run_lib
from absl import app, flags
from ml_collections.config_flags import config_flags
import os


FLAGS = flags.FLAGS

config_flags.DEFINE_config_file(
    'config', 'configs/training_ribodiffusion.py', 'Training configuration.', lock_config=True
)
flags.DEFINE_enum('mode', 'train', ['train', 'eval'],
                  'Running mode')
flags.DEFINE_integer('epochs', 10, 'Number of training epochs')
flags.DEFINE_float('learning_rate', 0.001, 'Learning rate for training')
flags.DEFINE_string('model_dir', 'checkpoints/', 'Directory to save trained model')
flags.DEFINE_boolean('deterministic', True, 'Set random seed for reproducibility.')

def main():
    create_model(FLAGS.config)
    if FLAGS.mode == 'inference':
        run_lib.vpsde_inference(FLAGS.config, FLAGS.save_folder, FLAGS.PDB_file)


if __name__ == '__name__':
    app.run(main)