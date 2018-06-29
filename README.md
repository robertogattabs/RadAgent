# RadAgent

## Preliminary Operations
First of all we have to load the library :)

```
library(RadAgent)
```

Then we have to load the CSV with the clinical outcome. The csv can also contains a set of clinical variable you may wish to add to the model. This is because, in general, in Radiomics investigation we are interested to find out the BEST model, which usually includes clinical and image features.
One note: the first line of the csv must contains the header

```
csv.Dino <- read.csv("~/original.csv.Dino.csv")
```
The first thing to do is to save the column with the clinical, dicotomic, outcome  (0/1), in an array and name the positions with the names of the folders of the corresponding DICOM Studies. In my case such names are numbers from 2 to 40 with some holes (I deleted some folders due to a not acceptable pixelSpacing)

```
arr.clinicalOutcome <- csv.Dino$Outcome
names(arr.clinicalOutcome) <- as.character( c(seq(2,6),seq(9,17),seq(19,24),seq(26,32),seq(34,40)) )
```

