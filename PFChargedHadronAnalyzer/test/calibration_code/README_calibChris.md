# For PFcalibration use PFCalibration/PFChargedHadronAnalyzer/test/calibration_code/calibChris.C  
Note: Please run calibChris.C code on ROOT version 5.xx (it might crash in ROOT version 6.xx)

Here are the following instructions :

1. You need to select only one _region_ from the lines at 39-42 in the calibChris code.

2. Then search "summary" in the code you will get list of commented plots. (you can find from lines 4368-4396)

    i)   There we have in total 12 drawGausFit functions for response wrt true energy. 
    
    ii)  There are 6 drawEtaDependence functions for response wrt abs(Eta).
    
    iii) There are 14 calibration coefficients plots for H barrel, H endcap, EH barrel & EH endcap.
    
   You need to uncomment lines to get the plot.
 
3. To run the code : 
```
root -l
root[0] .L calibChris.C+
root[1] calibChris()
```

For example :

If you want to look into the calibration coefficients for H barrel then you need to :
1. mention the _region_. ( char* _region_ = (char*)"barrel")
2. uncomment the H barrel calibration coefficients functions only.
Note: comment out all the other plots.    
3. Now, run and complile the code

If you want to look into the final corrected response wrt true energy for H barrel hadrons then you need to :
1. mention the _region_ . (char* _region_ = (char*)"barrel")
2. uncomment the drawGausFit(corrEtaBarrelHcal,responseCor,resolutionCor) function.
Note: comment out all the other plots.    
3. Now, run and complile the code and then you will get projections.root and output1.root along with the generated plots.

If you want to look into the final corrected response wrt eta for H hadrons then you need to :
1. mention the _region_ . (char* _region_ = (char*)"Full")
2. uncomment the drawEtaDependence(corrEtaDependenceH, responseEtaEtaH) function.
Note: comment out all the other plots.    
3. Now, run and complile the code and then you will get projections_eta.root and output2.root along with the generated plots.
