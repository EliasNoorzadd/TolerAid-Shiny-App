# DATA3888_KidneyA16

The following is a ReadMe file outlining the instructions on how to use the TolerAid shiny app. 

The folder named TolerAid should be downloaded from the github repo. Once downloaded, open the toleraid.Rproj project in Rstudio, and then open and run the app.R file.
These steps should launch the application window. In order to use the application correctly, the user should enter the relevant patient details as specified by the titles in the app, as well as a csv file of the patient's gene data in the correct format. The user will then select whether or not the patient is on immunotherapy. It is essential that the correct option is selected for if the patient is on immunotherapy, as this determines which model is used to predict the outcome. 
For Introduction purposes we have provided a folder which contains 4 sample csv files containing blood test data for 4 different patients. 
The following outline the names of the csv files and whether the patient associated with that file is on immunotherapy or not: 

non-TolerantFINAL.csv - not on immunotherapy,
tolerantNotOnImmunoFINAL.csv - not on immunotherapy,
tolerantOnImmunoFINAL.csv - on immunotherapy,
stableFINAL.csv - on immunotherapy,

By following the steps outlined above the user will understand how to use the app, and prevent the application from malfunctioning. In the event the application does crash or malfunction end the current R console from running, and run the app.R file again to continue as normal. 


