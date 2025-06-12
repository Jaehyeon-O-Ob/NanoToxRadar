# NanoToxRadar
## Introduction
NanoToxRadar is a multitarget nano-QSAR model with improved AD using a diverse multi-components nanoparticles (MC-NPs) dataset. This model was used to predict the cytotoxicity of MC-NP across 110 cell types.

Size-dependent electron-configuration fingerprint (SDEC FP) was employed to represent the structures of the MC-NPs[1], and one-hot encoded cell types were used to predict toxicities toward 110 cell lines. The CatBoost regression model achieved good performance $(R^{2}_{Test} = 0.877)$ and was deployed web-service. ([Link Text](https:www.kitox.re.kr/nanotoxradar))

If users want to know more information about it, please check the paper "NanoToxRadar: A multi-target nano-QSAR model for predicting the cytotoxicity of multicomponent nanoparticles" [2]

## Quick Usage


## Files
Each python code has explanation for how the code works.

Main code is prediction.py
if users want to use the model in users' local computer and follow the "Local"

## Local Usage
1. Frist Step - clone the repository
```bash
git clone https://github.com/Jaehyeon-O-Ob/NanoToxRadar.git
cd NanoToxRadar
```
2. Second Step - Create a virtual environment (Recommended with Miniconda)
```bash
conda create -n nanotoxradar python=3.10
```
3. Third Step - Install dependencies
```bash
pip install -r requirements.txt
```
4. Fourth Step - open 'prediction.py'
fill in the '# nano particle ready'
# nano particle ready
nanoparticle = {'Core':'CdSe',
               'Shell':'',
               'Doping': '',
               'Doping Rate(%)': '',
               'Coating': '',
               'Diameter(nm)': 500}
at this point, there are rules to enter each components.
- Core is needed.
- Shell is only filled with the list in "shell_volume_list.csv".
- If users input doping component, 'Doping Rage(%)' is mandatory within 1~100.
- The range of 'Diameter(nm)' covers from 1nm to 700nm. If nanomaterial diameter value is out of the range, the cytotoxicity prediction accuracy can be reduced as the model operates outside its trained domain.

5. Final Step - Users get "result_from_model.csv" file.

## References
[1] Shin et al., Use of Size-Dependent Electron Configuration Fingerprint to Develop General Prediction Models for Nanomaterials. NanoImpact 2021, 21, 100298.
[2] 
