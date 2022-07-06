#  Shell script to move files into proper location for docker image.  

#  This script assumes the docker image was created with a directory called output that can be mounted.  
#  In order to make this generic I use wild card characters and look for patterns.  
#  You can put multiple shell scripts and an overall shell script to run each individual experiment.  

cd /Prox_Labeling/data
mv *_Datafiles.zip /output/
mv *_Images.zip /output
