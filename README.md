<h1>RNAInvBench: Benchmarking RNA Inverse Folding with Secondary and Tertiary Structures</h1>
<p>RNAInvBench is a Python framework that provides environments, tasks, data, baseline algorithms, evaluation metrics and analysis tools for training/testing models for 2D and 3D RNA Inverse Design.
Notably within this benchmark, we split the overarching problem of RNA Inverse Design into three key tasks:

- **Secondary Pseudoknot-Free Inverse Design:** We aim to go from a 2D RNA structure, represented by dot-bracket notation, to a 1D RNA Sequence. We do not include the complex motif of pseudoknots, and most algorithms included here only include canonical base-pairings (A-U and G-C base-pairings only).
- **Secondary Pseudoknotted Inverse Design:** We aim to go from a 2D RNA structure that may or may not contain pseudoknots, to a 1D RNA Sequence. Notably the 2D RNA structures here are also represented by dot-bracket notation, but an extended version is used, where the type of bracket corresponds to the level of the pseudoknot. All pseudoknotted inverse design algorithms within this benchmark consider non-canonical pairings in their sequence design.
- **Tertiary Inverse Design:** We aim to go from a 3D RNA structure to a 1D RNA sequence. Key to this is to note that many 3D RNA structures are made up of a several chains of RNA, and thus they must be split into single chain RNAs before they can be properly fed into the inverse folding algorithms. The 3D RNA structure is represented by the PDB file, which is further converted into tensors for algorithmic representation.

  
The package contains the following algorithms:

- **Secondary Pseudoknot-Free Inverse Design:** 
  - antaRNA
  - aRNAque
  - GREED-RNA
  - IncaRNAtion
  - LEARNA
  - LibLEARNA
  - MCTS (Monte Carlo Tree Search)
  - Meta-LEARNA
  - Meta-LEARNA-Adapt
  - OmniGenome-GA
  - RNAInverse
  - SAMFEO
  - SentRNA
- **Secondary Pseudoknotted Inverse Design:** 
  - antaRNA
  - aRNAque
  - MCTS (Monte Carlo Tree Search)
- **Tertiary Inverse Design:** 
  - gRNAde
  - RDesign
  - RiboDiffusion

## Dataset Overview

### Training Data

| **Task**   | **Dataset**          | **Seq. & Struct. Count** | **Min Length** | **Max Length** | **Source**              |
|-----------|----------------------|-------------------------|---------------|---------------|-------------------------|
| PK-free  | Rfam-Learn-Train     | 65000                   | 50            | 450           | [Runge2019](https://arxiv.org/pdf/1812.11951) |
| PK-free  | TR0-PKfree           | 9815                    | 33            | 498           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| PK-inc   | TR0-PKinc            | 999                     | 37            | 496           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| PK-inc   | Pseudobase++         | 251                     | 21            | 137           | [Kleinkauf2015](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0815-6) |
| Tertiary | RNAsolo              | 4025                    | 11            | 4455          | [Joshi2023](https://arxiv.org/pdf/2305.14749) |

### Validation Data

| **Task**   | **Dataset**          | **Seq. & Struct. Count** | **Min Length** | **Max Length** | **Source**              |
|-----------|----------------------|-------------------------|---------------|---------------|-------------------------|
| PK-free  | VL0-PKfree           | 1193                    | 33            | 447           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| PK-inc   | VL0-PKinc            | 107                     | 53            | 497           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| Tertiary | RNAsolo              | 100                     | 11            | 714           | [Joshi2023](https://arxiv.org/pdf/2305.14749) |

### Testing Data

| **Task**   | **Dataset**            | **Seq. & Struct. Count** | **Min Length** | **Max Length** | **Source**              |
|-----------|------------------------|-------------------------|---------------|---------------|-------------------------|
| PK-free  | Eterna100-v1           | 100                     | 11            | 399           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| PK-free  | Eterna100-v2           | 100                     | 11            | 399           | [Koodli2021](https://www.biorxiv.org/content/10.1101/2021.08.26.457839v1) |
| PK-free  | Rfam-Learn-Test        | 100                     | 50            | 446           | [Runge2019](https://arxiv.org/pdf/1812.11951) |
| PK-free  | Rfam-Taneda            | 29                      | 54            | 451           | [Taneda2012](https://www.frontiersin.org/journals/genetics/articles/10.3389/fgene.2012.00036/full) |
| PK-free  | Rfam-Test              | 63                      | 35            | 273           | [Kleinkauf2015](https://academic.oup.com/bioinformatics/article/31/19/3114/210965) |
| PK-free  | RNA-Strand             | 50                      | 20            | 98            | [Kleinkauf2015](https://academic.oup.com/bioinformatics/article/31/19/3114/210965) |
| PK-free  | TS0-PKfree             | 1176                    | 27            | 499           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| PK-inc   | TS0-PKinc              | 129                     | 22            | 481           | [Singh2019](https://www.nature.com/articles/s41467-019-13395-9) |
| PK-inc   | Pseudobase++           | 251                     | 21            | 137           | [Kleinkauf2015](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-015-0815-6) |
| Tertiary | DAS-Split              | 98                      | 19            | 159           | [Joshi2023](https://arxiv.org/pdf/2305.14749) |
| Tertiary | structsim-v2-split     | 51                      | 21            | 373           | [Joshi2023](https://arxiv.org/pdf/2305.14749) |
| Tertiary | RNA-Puzzles            | 94                      | 6             | 256           | [Magnus2020RNAPuzzles](https://academic.oup.com/nar/article/48/2/576/5651330) |
| Tertiary | CASP15                 | 122                     | 43            | 720           | [Elofsson2023CASP15](https://www.sciencedirect.com/science/article/pii/S0959440X23000684) |

---

The package contains dockerfiles for each environment that is required to use the algorithms. Please make use of these when using the package!

We provide example code to run the algorithms, each example run uses the following naming convention: [algorithm_name]_run.py

</p>

<h2>Installation Instructions</h2>

<p>There are two key ways to run RNAInvBench. The recommended way is through docker, however you can also run it through creating your own environment through pyenv or conda.</p>

<h3>Docker Installation Instructions</h3>

<p>1. Install RNAInvBench source files.</p>
<p>2. Install <a href="https://docs.docker.com/get-docker/">docker</a> (if not done so already).</p>
<p>3. Navigate to the project directory.</p>
<p>4. Run <b>docker-compose up --build -d</b> to build the docker container.</p>
Once the docker-container has been fully built, you should see a new container in your docker application. You can then use the following commands to run the example main.py scripts provided:

<p>NOTE: RL represents the particular docker container to use, you will need to choose between rl, optim and antarna.</p>
<p>NOTE2: rnainvbench-main-rl-1 represents the name of the RL docker container, your name may differ - you should check the docker application for your own unique name.</p>
<p><b>docker exec -it $(docker-compose ps -q rl) rnainvbench-main-rl-1 python /app/rl_design_algos/main.py</b></p>

<h3>Conda Installation Instructions</h3>

<p>In each of the rl_design_algos, rna_design_algorithms and Tertiary_Design folders, an environment.yml file is provided. To run the algorithms within these folders, you need to install this environment.yml file.</p>
<p>Run <b>conda env create -f environment.yml</b> to install the environment.</p>
<p>To run the code afterwards, you simply need to use the main.py file provided within each of the folders.</p>
