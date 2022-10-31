This is the collection of code I used in the 2021 Year 3 Complexity & Networks Course at University. 
The first half of the course covered the Complexity of systems, and the concept of Self-Organised Critcality as demonstrated by the BTW and Oslo models.
The second half covered networks, considering models used to describe these networks and how they would develop over time.
Both halves of the course were assessed through a Project which required extensive use of python to create the data and analysis needed.
The code used for both projects is stored within this repository.

The Complexity Project looked at implementing the Oslo model, a simple model which explores Self-Organised Criticality. 
The model is iterative, beginnning with initialisation, before iterating through the drive and relaxation stages, which had to be coded within the model.
Aside from coding an implementation of the model, I also had to appropriately design it so that I could obtain information about the evolution of 
the system's variables over time, so that they could be analysed, and compared with different initialisation conditions. 
The functions and classes used are stored within Complexity_Module_Code.py, and they can be utilised to recreate the figures I used using the
ComplexityTesting.py script, which is designed for usage in the Spyder IDE, which uses #%% as a section break

The Networks Project considered three different variations of the Barabási-Albert (BA) model for random networks.
These models use different probabilities to determine whether a new node is connected to an old node within the network.
The Random Attachment variation determines this probability of this completely randomly.
The Preferential Attachment variation makes this probability proportional to the degree (no. of connections) of the old node.
Finally the Mixed Attachment variation has a pre-determined probability of using either of the aforementioned variations.
Utilising the Networkx library, I created a function which randomly generated a BA model with the ability to encode the variation of choice.
I also created a selection of functions to analyse the data I received to compare the relative performance of the different variations.
These functions are stored within Networks_Module_Code.py, and can also be used to recreate the data and figures I used using the NetworksTesting.py
script, which is also designed for use within the Spyder IDE.

Please Note also that the logbin code module was not created by myself, but provided as a tool to use for the projects, and thus is kept seperately here to
indicate its seperate origin.



