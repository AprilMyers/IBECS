# IBECS
IBECS is a pupil-analysis software for MATLAB. It takes in a ``PostCond`` file and calculates the following information: 
  * The area and dimensions of each pupil through the duration of every trial in the file
  * A graph of each pupil's dimensions contrasted against stimuli data
  * Trend behavior of the dimensions of the pupil as a response to stimuli

## Getting Started
These instructions will get you a copy of the project up and running on your local machine for development and testing purposes. 

## Prerequisites
This software is built for MATLAB and makes use of several built-in Add-Ons.

This software only works on data in the form of a ``PostCond`` file. In order to see the specifications of a ``PostCond`` file, see ``POSTCOND``. 

## Installing
In order to clone this GitHub repo, you can either use Command Line or download the files directly.

### For Windows (using Command Prompt):
```
to be added
```
### For Mac (using Terminal):
```
cd
cd Matlab/projects
git clone https://github.com/AprilMyers/IBECS.git
```

### Using direct download:
Navigate to ``https://github.com/AprilMyers/IBECS`` and click the green ``Clone or Download`` button. Download the ZIP file for the project. Open the ZIP file with a tool of the choice and put the resulting folder in your ``MATLAB/projects`` folder, or another location of your preference.  

## Running the Software
Open your downloaded project folder in MATLAB, and ensure that it is part of your file path. In MATLAB's command prompt, enter the command ``pupil``. You will be prompted to choose which parts of the data you want to have saved, rendered, and graphed. 

### Selecting a file
You will be prompted to choose what data the software will use. Use your computer's File Explorer/Finder to look for an appropriate ``PostCond`` file. 

### Cropping
After choosing the data you wish to be analyzed, you will be prompted to crop a frame of the data. The pupil selection requires two crops:
1. Select a bounding box around the pupil. A tighter bounding box around the pupil area will lead to better analysis. If you are satisfied with your cropping, follow the MATLAB command prompt to accept this crop. Otherwise, decline this crop and try again.
2. Select a bounding box for the parts of the pupil that are cut off. For example, if the top and bottom of the pupil are obfuscated by eyelids foobar.

## ``POSTCOND``
The ``PostCond`` data format is essentially a folder of TIF stacks and accompanying data for each trial, alongside a proper naming convention. The structure of your ``PostCond`` data should be as follows:
  * A TIF stack for each trial you wish to be analyzed (ENSURE that each trial is from the same session)
    * Each TIF stack should be named such that the string ``stimXXX`` is in it's name, where `XXX` represents the trial number. 
    * For example, the stack corresponding to trial 1 should contain the string `stim001`; the stack corresponding to trial 2 should contain the string `stim002`; and so on.
  * A ``.mat`` file that contains stimuli data for each trial. In particular
For each TIF stack, ensure the following conventions are met:
    * Update with specifications. 

## Authors
  * Ryan Natan - Postdoc Project Supervisor - UC Berkeley Systems of Neuroscience
  * April Myers - Graduate Student Advisor - UC Berkeley Vision Science    
  * Ethan Mehta - Undergraduate Research Apprentice - UC Berkeley MET
  * Varun Srivastava - Undergraduate Research Apprentice - UC Berkeley L&S