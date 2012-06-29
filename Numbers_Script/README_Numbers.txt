To produce the tex tables run ./Prediction_Producer.py 

Open up the script and change the path name of the files to the directory you output files are in

  .In the options at the top you specify the luminosity of the different samples
  .You can also specify the analysis 7/8TeV.
  . At the bottom of the script is the call to Number_Extractor in Bryn_Numbers.py which produces all of the tables. The options here are as follows
    . pass settings
    . sample dictionary
    . sample type ( One_btag/Two_btags/Zero_btags etc) You have to use the same naming convention that is currently thre
    . Apply Triggers True/False
    . Apply AlphaT in control samples (True/False) Zero btag is still run with/without so this will have to be run twice if you want to produce both tables
    . Calculation if specified will produce numbers via the formula method otherwise it will take the vanilla yields from the root files
    . Split Lumi/ Apply 3 diffrently lumis for the different samples / Apply single Lumi


#======================================================
To produce the closure tests run ./Prediction_Closure_Tests

Closure tests are computed within Jad_Compute.py, A dictionary is created and
passed to the Closure_Test function. 

A series of tests are defined within Jad_Compute. Addditional tests can be
added by adding a new entry i.e.

test_25 = {'MCS' : [], 'MCSE': [],'MCC': [], 'MCCE':[],'DC':[],'DS':[],'option' : -100.,'box' : "True", 'plot_title':"#mu + jets (#alpha_{T}>0.55) (0-btags) #rightarrow #mu#mu + jets (no #alpha_{T}) (0-btags)",'scale':None, 'reduce':"False", 'file_name':'Btag_mu_zero_to_dimuon_zero_alphaT_Cut','spread':"False"}

appending the new test to test_dicts and then including the parameters of the
closure test in ...

if self.file[self.entry]['AlphaT'] == '0.55' and self.file[self.entry]['Btag'] == 'Zero_btags' :
   self.Fill_Dictionary(test_25,Control = "Muon", Signal = "Muon",Not_Do = 'Signal') 
if self.file[self.entry]['AlphaT'] == '0.01' and self.file[self.entry]['Btag'] == 'Zero_btags' :
   self.Fill_Dictionary(test_25,Control = "DiMuon", Signal = "DiMuon",Not_Do = 'Control')  

#===================================================

To produce the Root Files for Teds limit code run
./Predition_Root_File_Producer.

A set of root files with the prefix RA1_Stats ....   are created. The trigger
corrections should be turned off when running the script as Ted implements
these corrections himself. 

