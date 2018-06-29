# RadAgent

## Preliminary Operations
First of all we have to load the library :)

```
library(RadAgent)
```

Then we have to load the CSV with the clinical outcome. The csv can also contains a set of clinical variable you may wish to add to the model. This is because, in general, in Radiomics investigation we are interested to find out the BEST model, which usually includes clinical and image features.
One note: the first line of the csv must contains the header

### Loading the clinical outcome
```
csv.Dino <- read.csv("~/original.csv.Dino.csv")
```
The first thing to do is to save the column with the clinical, dicotomic, outcome  (0/1), in an array and name the positions with the names of the folders of the corresponding DICOM Studies. In my case such names are numbers from 2 to 40 with some holes (I deleted some folders due to a not acceptable pixelSpacing)

```
arr.clinicalOutcome <- csv.Dino$Outcome
names(arr.clinicalOutcome) <- as.character( c(seq(2,6),seq(9,17),seq(19,24),seq(26,32),seq(34,40)) )
```

### Loading other clinical covariates
This step is optional: you can train a Radiomics model based only on features acquired from the DICOM images.
However, as example, now I load other clinical covariates, such as Age, Sex, status of the pelivic lymphnodes, cT and cN. I build a new matrix called 'other.covariates.matrix':

```
# assign the AGE, present in the csv (numeric)
arr.eta <- csv.Dino$eta
# assign teh SEX, present in csv (0,1)
arr.sesso <- csv.Dino$sesso
# assign the status of pelvic lynfonode, present in the csv (0,1)
arr.STATO.LINFONODI.INGUINALI <- csv.Dino$STATO.LINFONODI.INGUINALI
# assigne the clinical T stage, present in the CSV (categorical 1,2,3,4)
arr.xT <- csv.Dino$xT
# assigne the clinical N stage, present in the CSV (categorical 1,2,3)
arr.xN <- csv.Dino$N
# build the matrix with the clinical variables to include in the model
other.covariates.matrix <- cbind( "eta"=arr.eta, "sesso"=arr.sesso, 
                                  "STATO.LINFONODI.INGUINALI"=arr.STATO.LINFONODI.INGUINALI, 
                                  "xT"=arr.xT, "xN"=arr.xN )                                 
```

### RADIOMICS agent
RadAgent extends the moddicom package [1], providing the entire pipeline for Radiomics investigation. This means that RadAgent:
* extract the image features on the found foldes, given and entry point on the filesystem. The features are compliant with the Image biomarker standardisation initiative [2].
* given and entry point on the filesystem, it loads the images from the found foldes. 
* applies the LoGs filtering, according to the given sigmas;
* for each filtered image, it extracts the features, compliant with the Image biomarker standardisation initiative [2].
* for each feature, changing the sigmas in the LoG filtering, it calculate makes a Mann Withney test with the clinical outcome. Then it takes the sigmas able to give the highest performances to that feature. In case of high variation of the p.value for small variations of the sigma, it can discharge the feature (doubt of overfitting)
* it check the couple of features and in case of high correlation (pearson test) it discharge the feature with the lowest performance in predicting the results;
* it begins to build a model with a feedforeward feature selection strategy, until the number of requested feature in the model is reached.


---
references:

[1] Dinapoli N, Alitto AR, Vallati M, Gatta R, Autorino R, Boldrini L, Damiani A, Valentini V. Moddicom: a complete and easily accessible library for prognostic evaluations relying on image features. Conf Proc IEEE Eng Med Biol Soc. 2015 Aug;2015:771-4. doi: 10.1109/EMBC.2015.7318476.

[2] Zwanenburg A, Leger S, Vallières M, Löck S. Image biomarker standardisation initiative, eprint arXiv:1612.07003. Available (29/06/2018) at: https://arxiv.org/abs/1612.07003 


---





