#!cafe
#version
#date

load -i taxon.tab -l logfile.txt -p 0.05
tree (((((Craseoa_lathetica:17,Abylopsis_tetragona:17):7,(Agalma_elegans:14,Nanomia_bijuga:14):10):7,Physalia_physalis:32):19,Hydra_magnipapillata:51):48,Nematostella_vectensis:100);
lambdamu -s -t
report resultfile