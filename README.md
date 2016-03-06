# BCML to SBML converter

The Biological Connection Markup Language (BCML) format, as first described by Beltrame, Calura, *et al.* (2011), has been successfully used in the past years in order to describe multiple regulatory networks and pathways, especially Dendritic Cells-related intracellular signalling pathways. These pathways have been collected into a set of maps available in the [DC-Atlas resource](http://compbiotoolbox.fmach.it/DCATLAS.php).

In order to visualise this data using softwares similar to CellDesigner, we developed a simple converter from BCMl to SBML (version 2, level 4).


## Workflow

**Step 1** : Convert from BCMl to SBML

Read a BCML file and convert it into SBML (version 2, level 4) using the bcml_to_sbml.py script. This step uses SBO terms in order not to lose reaction annotations.

  $ python bcml_to_sbml.py TLR9.xml

**Step 2** : Open newly created files in CellDesigner

Open newly created files in CellDesigner and save them under a different name. CellDesigner will add its extended content to the file. This extended content can then be adjusted in the next step to take advantage of all CellDesigner functionalities and visual representations.

**Step 3** : Improve visualisation

Improve XML files generated in Step 3 by CellDesigner, using improve_cd_file.py.

  $ python improve_cd_file.py TLR9_cd.xml

**Step 4** : Manual adjustments

Once the conversion has been made, users can reorganise to their liking the newly created maps using CellDesigner.

