The analysis code should run entirely out of the box.

To change pileup weights go into batchgoldensinglemu.py and search out PU_2012. To add and new json go to the last line of the script and add in the path to the json file.

Simply run ./Run_All_MC ./Run_All_Data

Go to the output directory and then run ./Make_All_Muon/Had/Photon once the jobs are completed.

Next transfer the output root files to somewhere where you can run the plotting script/Formula method code
