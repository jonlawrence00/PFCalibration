## ./run.sh <path of folderr containing all output files>
path=$1

#### To get validation and training loss wrt epoch
#./run.sh output/target_trueE/Trained_2_350_eta3_epoch100/


python3.10 scripts/lossplot.py ${path}
#### To get output.root file having 1D energy responses, response & resolution wrt E(true) inclusive in eta
#python3.10 full_binned.py ${path} output.root ${path}pickle_test/
mkdir ${path}/pdf
rm -rf ${path}/pdf
mkdir ${path}/pdf
cp -r res index.* ${path}/pdf/.

for i in EH_barrel EH_ec_in EH_ec_out H_barrel H_ec_in H_ec_out
do
#    rm ${path}/pdf/${i}/*
    mkdir ${path}/pdf/${i}
    cp -r res index.* ${path}/pdf/${i}/.
done



root -l -q 'scripts/genrateplot_responsevsE_separate.C("EH (0 < |#eta| <1.55)","barrel_corrEtaBarrelEcalHcal","'${path}'","output_eta2pt75_EH","barrel")'
root -l -q  'scripts/genrateplot_responsevsE_separate.C("EH hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","ec_in")'
root -l -q  'scripts/genrateplot_responsevsE_separate.C("EH hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","ec_out")'

root -l -q  'scripts/genrateplot_responsevsE_separate.C("H (0 < |#eta| <1.55)","barrel_corrEtaBarrelHcal","'${path}'","output_eta2pt75_H","barrel")'
root -l -q  'scripts/genrateplot_responsevsE_separate.C("H hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","ec_in")'
root -l -q  'scripts/genrateplot_responsevsE_separate.C("H hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","ec_out")'

root -l -q 'scripts/genrateplot_resovsE_separate.C("EH (0 < |#eta| <1.55)","barrel_corrEtaBarrelEcalHcal","'${path}'","output_eta2pt75_EH","barrel")'
root -l -q  'scripts/genrateplot_resovsE_separate.C("EH hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","ec_in")'
root -l -q  'scripts/genrateplot_resovsE_separate.C("EH hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","ec_out")'
root -l -q  'scripts/genrateplot_resovsE_separate.C("H (0 < |#eta| <1.55)","barrel_corrEtaBarrelHcal","'${path}'","output_eta2pt75_H","barrel")'
root -l -q  'scripts/genrateplot_resovsE_separate.C("H hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","ec_in")'
root -l -q  'scripts/genrateplot_resovsE_separate.C("H hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","ec_out")'
root -l -q 'scripts/genrateplot_responsevsEta.C("EH hadrons","corrEtaDependenceEH_wrtEta","'${path}'","output_eta2pt75_EH_eta","EH_all")'
root -l -q 'scripts/genrateplot_responsevsEta.C("H hadrons","corrEtaDependenceH_wrtEta","'${path}'","output_eta2pt75_H_eta","H_all")'

root -l -q  'scripts/genrateplot_1dresponse_eta.C("EH hadrons","corrEtaDependenceEH","'${path}'","output_eta2pt75_EH_eta","trans_EH")'
root -l -q  'scripts/genrateplot_1dresponse_eta.C("H hadrons","corrEtaDependenceH","'${path}'","output_eta2pt75_H_eta","trans_H")'


#root -l -q 'scripts/genrateplot_responsevsE_UL2018wrtrun3.C("EH (0 < |#eta| <1.55)","barrel_corrEtaBarrelEcalHcal","'${path}'","output_eta2pt75_EH","barrel")'

#root -l -q  'scripts/genrateplot_responsevsE_UL2018wrtrun3.C("EH hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","ec_in")'
#root -l -q  'scripts/genrateplot_responsevsE_UL2018wrtrun3.C("EH hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","ec_out")'

#root -l -q  'scripts/genrateplot_responsevsE_UL2018wrtrun3.C("H (0 < |#eta| <1.55)","barrel_corrEtaBarrelHcal","'${path}'","output_eta2pt75_H","barrel")'
#root -l -q  'scripts/genrateplot_responsevsE_UL2018wrtrun3.C("H hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","ec_in")'
#root -l -q  'scripts/genrateplot_responsevsE_UL2018wrtrun3.C("H hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","ec_out")'



root -l -q 'scripts/genrateplot_1dresponsevsE_UL2018wrtrun3.C("EH (0 < |#eta| <1.55)","barrel_corrEtaBarrelEcalHcal","'${path}'","Output_eta2pt75_EH","EH_barrel",20,24)'
root -l -q  'scripts/genrateplot_1dresponsevsE_UL2018wrtrun3.C("EH hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","EH_ec_in")'
root -l -q  'scripts/genrateplot_1dresponsevsE_UL2018wrtrun3.C("EH hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapEcalHcal","'${path}'","output_eta2pt75_EH","EH_ec_out")'
root -l -q  'scripts/genrateplot_1dresponsevsE_UL2018wrtrun3.C("H (0 < |#eta| <1.55)","barrel_corrEtaBarrelHcal","'${path}'","output_eta2pt75_H","H_barrel")'
root -l -q  'scripts/genrateplot_1dresponsevsE_UL2018wrtrun3.C("H hadrons (1.55 < |#eta| <2.5)","EC_within_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","H_ec_in")'
root -l -q  'scripts/genrateplot_1dresponsevsE_UL2018wrtrun3.C("H hadrons (2.5 < |#eta| <2.75)","EC_outside_tracker_corrEtaEndcapHcal","'${path}'","output_eta2pt75_H","H_ec_out")'

####
