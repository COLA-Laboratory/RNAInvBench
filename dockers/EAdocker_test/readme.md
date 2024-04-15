## note
antaRNA is originally write in Python2.x, I used `2to3` tool to convert them into python3 version code first. So files in `antarna` is slightly different from original code of the paper.

## how to use
```shell
cd RNADesignPipeline/dockers/EAdocker_test  # work dir on your local machine
docker build -t rnabench:0.5.0 .  # build image. '0.5.0' is my personal version flag, you can change it, first build will take some time
docker run -it --rm rnabench:0.5.0  # run container from image. image and version flag should be identical to the build command. --rm means delete this container after closing it, you can also save this container by removing this option
# in container
cd /rnabench  # work dir I created in image, and you can see you are in a conda environment named rnabench
# test antaRNA using file: antaRNA_v109.py(provided by antaRNA group)
python antaRNA_v109.py -h
# you should see some help informations :)
```


