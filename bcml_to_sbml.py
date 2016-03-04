#!/usr/bin/python

# Read a BCML file and create an SBML.
# 
# Example of a command line : 
# python bcml_to_sbml.py ../DC-ATLAS/BCML/TLR9.xml
# 

# General
import sys
import os.path
import re

# For BCML
import xml.etree.ElementTree as ET

# For SBML
from libsbml import *

def check(value, message):
	"""If 'value' is None, prints an error message constructed using
	'message' and then exits with status code 1.  If 'value' is an integer,
	it assumes it is a libSBML return status code.  If the code value is
	LIBSBML_OPERATION_SUCCESS, returns without further action; if it is not,
	prints an error message constructed using 'message' along with text from
	libSBML explaining the meaning of the code, and exits with status code 1.
	"""
	if value == None:
		raise SystemExit('LibSBML returned a null value trying to ' + message + '.')
	elif type(value) is int:
		if value == LIBSBML_OPERATION_SUCCESS:
			return
		else:
			err_msg = 'Error encountered trying to ' + message + '.' \
				+ 'LibSBML returned error code ' + str(value) + ': "' \
				+ OperationReturnValue_toString(value).strip() + '"'
			raise SystemExit(err_msg)
	else:
		return


def addMacroMolecule(bcmlMacromol, sbmlModel, sbmlCompartmentId):
	
	sbmlSpecies = sbmlModel.createSpecies()
	#print(bcmlMacromol.attrib.get('ID'))
	
	# Id
	check( sbmlSpecies.setId(idfy(str(bcmlMacromol.attrib.get('ID'))) ),				"Set ID")
	# Name
	label = str(bcmlMacromol.attrib.get('label'))
	if label is None or label == '' or label == 'None':
		label = str(bcmlMacromol.attrib.get('cloneref'))
	check( sbmlSpecies.setName(label),					"Set name")
	# Compartment
	check( sbmlSpecies.setCompartment(sbmlCompartmentId),								"Set compartment")
	# SpeciesType
	#check( sbmlSpecies.setSpeciesType("protein"),										"Set compartment")
	
	# Add note
	check( sbmlSpecies.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\n"+extractNotes(bcmlMacromol)+"</p>"), 							"Add notes")

def addNucleicAcidFeature(bcmlNAfeature, sbmlModel, sbmlCompartmentId):
	
	sbmlSpecies = sbmlModel.createSpecies()
	#print(bcmlNAfeature.attrib.get('ID'))
	
	# Id
	check( sbmlSpecies.setId(idfy(str(bcmlNAfeature.attrib.get('ID'))) ),				"Set ID")
	# Name
	check( sbmlSpecies.setName(str(bcmlNAfeature.attrib.get('label'))),					"Set name")
	# Compartment
	check( sbmlSpecies.setCompartment(sbmlCompartmentId),								"Set compartment")
	# SpeciesType
	#if bcmlNAfeature.find('UnitOfInformation') is not None:
	#	naFeatureLabel =  bcmlNAfeature.find('UnitOfInformation').attrib.get('label')
	#	if naFeatureLabel == 'mRNA':
	#		check( sbmlSpecies.setSpeciesType("RNA"),									"Set type")
	#	if naFeatureLabel == 'gene':
	#		check( sbmlSpecies.setSpeciesType("gene"),									"Set type")
	
	# Add note
	check( sbmlSpecies.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\n"+extractNotes(bcmlNAfeature)+"</p>"), 							"Add notes")

def addSimpleChemical(bcmlSimpleChem, sbmlModel, sbmlCompartmentId):
	
	sbmlSpecies = sbmlModel.createSpecies()
	#print(bcmlSimpleChem.attrib.get('ID'))
	
	# Id
	check( sbmlSpecies.setId(idfy(str(bcmlSimpleChem.attrib.get('ID')))),					"Set ID")
	# Name
	check( sbmlSpecies.setName(str(bcmlSimpleChem.attrib.get('label'))),					"Set name")
	# Compartment
	check( sbmlSpecies.setCompartment(sbmlCompartmentId),									"Set compartment")
	# SpeciesType
	#check( sbmlSpecies.setSpeciesType("simple"),											"Set type")
	# SBO:0000247 - simple chemical
	check( sbmlSpecies.setSBOTerm(247),															"Set SBO term")	
	
	# Add note
	check( sbmlSpecies.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\n"+extractNotes(bcmlSimpleChem)+"</p>"), 								"Add notes")

def addComplex(bcmlComplex, sbmlModel, sbmlCompartmentId):
	
	sbmlSpecies = sbmlModel.createSpecies()
	#print(idfy(str(bcmlComplex.attrib.get('ID'))))
	
	# Id
	check( sbmlSpecies.setId(idfy(str(bcmlComplex.attrib.get('ID')))),						"Set ID")
	# Name
	check( sbmlSpecies.setName(str(bcmlComplex.attrib.get('ID'))+":"),						"Set name")
	# Compartment
	check( sbmlSpecies.setCompartment(sbmlCompartmentId),									"Set compartment")
	# SpeciesType
	#check( sbmlSpecies.setSpeciesType("complex"),											"Set type")
	
	# Add notes
	supplementaryNotes = ""
	
	# Add cardinality information as note (if present)
	if bcmlComplex.attrib.get('cardinality') is not None:
		supplementaryNotes += "Complex:Cardinality:"+str(bcmlComplex.attrib.get('cardinality'))+"\n"
	
	# Add supplementary macromolecules as species
	countNewMacroMolecules = 0
	countAll = 0
	for bcmlMacromol in bcmlComplex.findall('Macromolecule'):
		countAll+=1
		if bcmlMacromol.attrib.get('cloneref') is None:
			addMacroMolecule(bcmlMacromol, sbmlModel, sbmlCompartmentId)
			countNewMacroMolecules+=1
			supplementaryNotes += "Complex:MacroMolecule:"+idfy(str(bcmlMacromol.attrib.get('ID')))+"\n"
		else:
			supplementaryNotes += "Complex:MacroMolecule:"+idfy(str(bcmlMacromol.attrib.get('cloneref')))+"\n"
	
	# Add supplementary complex as species (yes, there can be complexes in complexes)
	for bcmlSubComplex in bcmlComplex.findall('Complex'):
		countAll+=1
		if bcmlSubComplex.attrib.get('cloneref') is None:
			addComplex(bcmlSubComplex, sbmlModel, sbmlCompartmentId)
			countNewMacroMolecules+=1
			supplementaryNotes += "Complex:Complex:"+idfy(str(bcmlSubComplex.attrib.get('ID')))+"\n"
		else:
			supplementaryNotes += "Complex:Complex:"+idfy(str(bcmlSubComplex.attrib.get('cloneref')))+"\n"
	
	# Can also contain SimpleChemicals
	for bcmlSimpleChem in bcmlComplex.findall('SimpleChemical'):
		countAll+=1
		if bcmlSimpleChem.attrib.get('cloneref') is None:
			addSimpleChemical(bcmlSimpleChem, sbmlModel, sbmlCompartmentId)
			countNewMacroMolecules+=1
			supplementaryNotes += "Complex:SimpleChemical:"+idfy(str(bcmlSimpleChem.attrib.get('ID')))+"\n"
		else:
			supplementaryNotes += "Complex:SimpleChemical:"+idfy(str(bcmlMacromol.attrib.get('cloneref')))+"\n"
	
	# Create a new association reaction if all macromolecules inside the complex are newly defined
	if countAll==countNewMacroMolecules:
		sbmlReaction = sbmlModel.createReaction()
		# Add notes saying to report to the complex notes
		check( sbmlReaction.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\nReaction:Complex building\nSee produced Complex for details.</p>"),	"Add notes")
		
		if bcmlComplex.attrib.get('type') == 'And':
			# ID
			check( sbmlReaction.setId("ra"+idfy(str(bcmlComplex.attrib.get('ID')))),											"Set ID")
			# Is reversible
			check( sbmlReaction.setReversible(False),														"Set reversible")
			# SBO term: and
			check( sbmlReaction.setSBOTerm(173),															"Set SBO term")	
			#print(sbmlReaction.getId())
		
		elif bcmlComplex.attrib.get('type') == 'Or':
			# ID
			check( sbmlReaction.setId("ro"+idfy(str(bcmlComplex.attrib.get('ID')))),											"Set ID")
			# Is reversible
			check( sbmlReaction.setReversible(False),														"Set reversible")
			# SBO term: or
			check( sbmlReaction.setSBOTerm(174),															"Set SBO term")	
			#print(sbmlReaction.getId())
			
		else:
			# ID
			check( sbmlReaction.setId("ru"+idfy(str(bcmlComplex.attrib.get('ID')))),											"Set ID")
			# Is reversible
			check( sbmlReaction.setReversible(False),														"Set reversible")
			# SBO term: logical combination (logic unknown or not specified)
			check( sbmlReaction.setSBOTerm(237),															"Set SBO term")	
			#print(sbmlReaction.getId())

		# Add reactants
		for bcmlReactant in bcmlComplex.findall('Macromolecule'):
			sbmlReactant = sbmlReaction.createReactant()
			#print("- "+bcmlReactant.attrib.get('label')+" "+idfy(bcmlReactant.attrib.get('ID')))
			check( sbmlReactant.setSpecies(idfy(bcmlReactant.attrib.get('ID'))),					"Set Reference")
		
		# Add product
		sbmlProduct = sbmlReaction.createProduct()
		#print("+ "+idfy(bcmlComplex.attrib.get('ID')))
		check( sbmlProduct.setSpecies(idfy(bcmlComplex.attrib.get('ID'))),						"Set Reference")
	
	# Include logic information in notes
	if bcmlComplex.attrib.get('type') == 'And':
		supplementaryNotes += "Complex:Logic:And\n"
	elif bcmlComplex.attrib.get('type') == 'Or':
		supplementaryNotes += "Complex:Logic:Or\n"
	else:
		supplementaryNotes += "Complex:Logic:?\n"

	# Add note
	check( sbmlSpecies.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\n"+supplementaryNotes+"\n"+extractNotes(bcmlComplex)+"</p>"),	"Add notes")

def addProcess(bcmlReaction, sbmlModel, sbmlCompartmentId, sbmlReactionNb, andDict, orDict, bcmlRoot):
    
	bcmlReactant = bcmlReaction.find('Consumption')
	bcmlProduct = bcmlReaction.find('Production')
	
	supplementaryNotes = ""
	
	# Create reaction
	sbmlReaction = sbmlModel.createReaction()
	# Is reversible
	check( sbmlReaction.setReversible(False),													"Set reversible")

	# Checking if we have a case of transcription Source -> mRNA (with gene and TF referenced through a AndNode)
	if str(bcmlProduct.attrib.get('refNode')).startswith("mRNA") and re.match("^[Ss][0-9]{1,2}$", str(bcmlReactant.attrib.get('refNode')), flags=0) is not None:
		
		supplementaryNotes += "Reaction:Transcription\n"
		
		# ID
		check( sbmlReaction.setId("tra"+str(sbmlReactionNb)),										"Set ID")
		#print(sbmlReaction.getId())
		
		# SBO:0000183 - transcription
		check( sbmlReaction.setSBOTerm(183),															"Set SBO term")	
		
		bcmlStimulation = bcmlReaction.find('NecessaryStimulation')
		bcmlStimulationRefNode = bcmlStimulation.attrib.get('refNode')
		
		if bcmlStimulationRefNode in andDict:
			
			bcmlAndNode = bcmlRoot.find(".//AndNode[@ID='"+bcmlStimulationRefNode+"']")
			for bcmlLog in bcmlAndNode.findall('Logic'):
				if bcmlLog.attrib.get('refNode').startswith('gene'):
					supplementaryNotes += addReactant(bcmlLog, sbmlReaction, orDict, bcmlRoot)
				else:
					supplementaryNotes += addModulation(bcmlLog, sbmlReaction, orDict, bcmlRoot)

			# RNA
			supplementaryNotes += addProduct(bcmlProduct, sbmlReaction, orDict, bcmlRoot)

		else:
			
			# Gene=>reactant, or else=>modulation?
			if bcmlStimulationRefNode.startswith('gene'):
				supplementaryNotes += addReactant(bcmlStimulation, sbmlReaction, orDict, bcmlRoot)
			else:
				supplementaryNotes += addModulation(bcmlStimulation, sbmlReaction, orDict, bcmlRoot)
			# RNA
			supplementaryNotes += addProduct(bcmlProduct, sbmlReaction, orDict, bcmlRoot)			
	else:
		
		supplementaryNotes += "Reaction:Generic\n"
		
		# ID
		check( sbmlReaction.setId("re"+str(sbmlReactionNb)),										"Set ID")
		#print(sbmlReaction.getId())

		# Add reactants
		for bcmlReactant in bcmlReaction.findall('Consumption'):
			supplementaryNotes += addReactant(bcmlReactant, sbmlReaction, orDict, bcmlRoot)
		
		# Add product
		for bcmlProduct in bcmlReaction.findall('Production'):	
			supplementaryNotes += addProduct(bcmlProduct, sbmlReaction, orDict, bcmlRoot)
		
		# Add necessary stimulation
		# SBO:0000461 - essential activator
		for bcmlNecStim in bcmlReaction.findall('NecessaryStimulation'):	
			supplementaryNotes += addNecessaryStimulation(bcmlNecStim, sbmlReaction, orDict, bcmlRoot)

	
	# Add modulation
	# SBO:0000462 - non essential stimulator
	for bcmlModifier in bcmlReaction.findall('Modulation'):	
		supplementaryNotes += addModulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add inhibitor
	# SBO:0000020 - inhibitor
	for bcmlModifier in bcmlReaction.findall('Inhibition'):	
		supplementaryNotes += addInhibition(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add catalysis
	# SBO:0000013 - catalyst
	for bcmlModifier in bcmlReaction.findall('Catalysis'):	
		supplementaryNotes += addCatalysis(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add stimulation
	# SBO:0000459 - stimulator
	for bcmlModifier in bcmlReaction.findall('Stimulation'):	
		supplementaryNotes += addStimulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add note
	check( sbmlReaction.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\n"+supplementaryNotes+"</p>"),	"Add notes")

def addReactant(bcmlReactant, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlReactant.attrib.get('refNode')
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addReactant(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
	
	# Reactant should not be a source
	if re.match("^[Ss][0-9]{1,2}$", str(speciesRef), flags=0) is None:
		sbmlReactant = sbmlReaction.createReactant()
		check( sbmlReactant.setSpecies(idfy(speciesRef)),											"Set Reference")
		#print("- "+speciesRef+" "+idfy(speciesRef))
		return "Reactant:"+idfy(speciesRef)+"\n"
		
	return ""
	
def addProduct(bcmlProduct, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlProduct.attrib.get('refNode')
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addProduct(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
		
	# Product should not be a sink
	if re.match("^[Ss][0-9]{1,2}$", str(speciesRef), flags=0) is None:
		sbmlProduct = sbmlReaction.createProduct()
		check( sbmlProduct.setSpecies(idfy(speciesRef)),											"Set Reference")
		#print("+ "+speciesRef+" "+idfy(speciesRef))
		return "Product:"+idfy(speciesRef)+"\n"
	
	return ""

def addModulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlModifier.attrib.get('refNode')
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addModulation(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
		
	sbmlModifier = sbmlReaction.createModifier()
	check( sbmlModifier.setSpecies(idfy(speciesRef)),												"Set Reference")
	# SBO:0000462 - non essential stimulator
	check( sbmlModifier.setSBOTerm(462),															"Set SBO term")	
	
	#print("o "+speciesRef+" "+idfy(speciesRef))
	return "Modulation:"+idfy(speciesRef)+"\n"

def addInhibition(bcmlModifier, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlModifier.attrib.get('refNode')
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addInhibition(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
		
	sbmlModifier = sbmlReaction.createModifier()
	check( sbmlModifier.setSpecies(idfy(speciesRef)),												"Set Reference")
	# SBO:0000020 - inhibitor
	check( sbmlModifier.setSBOTerm(20),																"Set SBO term")	
	
	#print("x "+speciesRef+" "+idfy(speciesRef))
	return "Inhibition:"+idfy(speciesRef)+"\n"

def addCatalysis(bcmlModifier, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlModifier.attrib.get('refNode')
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addCatalysis(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
	
	sbmlModifier = sbmlReaction.createModifier()
	check( sbmlModifier.setSpecies(idfy(speciesRef)),												"Set Reference")
	# SBO:0000013 - catalyst
	check( sbmlModifier.setSBOTerm(13),																"Set SBO term")	
	
	#print("c "+speciesRef+" "+idfy(speciesRef))
	return "Catalysis:"+idfy(speciesRef)+"\n"

def addNecessaryStimulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlModifier.attrib.get('refNode')
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addNecessaryStimulation(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
		
	sbmlModifier = sbmlReaction.createModifier()
	check( sbmlModifier.setSpecies(idfy(speciesRef)),												"Set Reference")
	# SBO:0000461 - essential activator
	check( sbmlModifier.setSBOTerm(461),															"Set SBO term")	
	
	#print("S "+speciesRef+" "+idfy(speciesRef))
	return "NecessaryStimulation:"+idfy(speciesRef)+"\n"

def addStimulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot):
	speciesRef = bcmlModifier.attrib.get('refNode')
	if speciesRef is None or speciesRef == '':
		speciesRef = bcmlModifier.text
	if speciesRef is None:
		raise SystemExit('Stimulation without a reference! ' + speciesRef + '.')
	
	if speciesRef in orDict:
		bcmlOrNode = bcmlRoot.find(".//OrNode[@ID='"+speciesRef+"']")
		notes = ""
		for bcmlLog in bcmlOrNode.findall('Logic'):
			notes += addStimulation(bcmlLog, sbmlReaction, orDict, bcmlRoot)
		return notes
		
	sbmlModifier = sbmlReaction.createModifier()
	check( sbmlModifier.setSpecies(idfy(speciesRef)),												"Set Reference")
	# SBO:0000459 - stimulator
	check( sbmlModifier.setSBOTerm(459),															"Set SBO term")	
	
	#print("s "+speciesRef+" "+idfy(speciesRef))
	return "Stimulation:"+idfy(speciesRef)+"\n"

def addReaction(bcmlReaction, sbmlModel, sbmlCompartmentId, sbmlReactionNb, orDict, bcmlRoot):
	"""Association or dissociation"""
	
	sbmlReaction = sbmlModel.createReaction()

	# ID
	check( sbmlReaction.setId("re"+str(sbmlReactionNb)),											"Set ID")
	# Compartment: can't set compartment for SMBL2.4
	#check( sbmlReaction.setCompartment(sbmlCompartmentId),											"Set compartment")
	# Is reversible
	check( sbmlReaction.setReversible(False),														"Set reversible")
	#print(sbmlReaction.getId())

	# Add notes
	supplementaryNotes = ""
	
	# Keep reaction type in notes
	supplementaryNotes += "Reaction:"+bcmlReaction.tag+"\n"
		
	# Add reactants
	for bcmlReactant in bcmlReaction.findall('Consumption'):
		supplementaryNotes += addReactant(bcmlReactant, sbmlReaction, orDict, bcmlRoot)

	# Add product
	for bcmlProduct in bcmlReaction.findall('Production'):	
		supplementaryNotes += addProduct(bcmlProduct, sbmlReaction, orDict, bcmlRoot)
	
	# Add modulation
	# # SBO:0000462 - non essential stimulator
	for bcmlModifier in bcmlReaction.findall('Modulation'):	
		supplementaryNotes += addModulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
		
	# Add inhibitor
	# SBO:0000020 - inhibitor
	for bcmlModifier in bcmlReaction.findall('Inhibition'):	
		supplementaryNotes += addInhibition(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
		
	# Add catalysis
	# SBO:0000013 - catalyst
	for bcmlModifier in bcmlReaction.findall('Catalysis'):	
		supplementaryNotes += addCatalysis(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add necessary stimulation
	# SBO:0000461 - essential activator
	for bcmlModifier in bcmlReaction.findall('NecessaryStimulation'):	
		supplementaryNotes += addNecessaryStimulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add stimulation
	# SBO:0000459 - stimulator
	for bcmlModifier in bcmlReaction.findall('Stimulation'):	
		supplementaryNotes += addStimulation(bcmlModifier, sbmlReaction, orDict, bcmlRoot)
	
	# Add note
	check( sbmlReaction.setNotes("<p xmlns='http://www.w3.org/1999/xhtml'>\n"+supplementaryNotes+"</p>"),	"Add notes")


def extractNotes(bcmlElement):
	hasNotes = False
	sbmlNotes = ""
	
	# All children of Finding become note lines
	if bcmlElement.find('Finding') is not None:
		hasNotes = True
		# Loop through all Finding children
		for bcmlFinding in bcmlElement.find('Finding'):
			if bcmlFinding.text is not None:
				sbmlNotes += bcmlFinding.tag+":"+bcmlFinding.text.strip()+"\n"
	
	# MacroModule become note line
	# e.g. MacroModule:ReceptorSensing
	if bcmlElement.find('MacroModule') is not None:
		hasNotes = True
		# Loop through all MacroModules
		for bcmlMacroModule in bcmlElement.findall('MacroModule'):
			if bcmlMacroModule.text is not None:
				# MODULE:xxx ?
				sbmlNotes += bcmlMacroModule.tag+":"+bcmlMacroModule.text.strip()+"\n"
	
	# StateVariable become note line
	# e.g. MacroModule:ReceptorSensing
	if bcmlElement.find('StateVariable') is not None:
		hasNotes = True
		# Loop through all StateVariable
		for bcmlStateVariable in bcmlElement.findall('StateVariable'):
			labelTxt = bcmlStateVariable.attrib.get('label')
			if labelTxt is not None and labelTxt is not '' and labelTxt is not ' ':
				# label="inactive", label="P@507"
				sbmlNotes += bcmlStateVariable.tag+":"+bcmlStateVariable.attrib.get('label')+"\n"
	
	# UnitOfInformation become note line
	# e.g. MacroModule:ReceptorSensing
	if bcmlElement.find('UnitOfInformation') is not None:
		hasNotes = True
		# Loop through all StateVariable
		for bcmlUnitOfInfo in bcmlElement.findall('UnitOfInformation'):
			# label="open" or "close"
			# mRNA, gene
			if bcmlUnitOfInfo.attrib.get('label') is not None:
				sbmlNotes += bcmlUnitOfInfo.tag+":label:"+bcmlUnitOfInfo.attrib.get('label')+"\n"
			if bcmlUnitOfInfo.attrib.get('prefix') is not None:
				# <UnitOfInformation prefix="mt" term="psac"/>
				sbmlNotes += bcmlUnitOfInfo.tag+":prefix:"+bcmlUnitOfInfo.attrib.get('prefix')+"\n"
			if bcmlUnitOfInfo.attrib.get('term') is not None:
				# <UnitOfInformation prefix="mt" term="psac"/>
				sbmlNotes += bcmlUnitOfInfo.tag+":term:"+bcmlUnitOfInfo.attrib.get('term')+"\n"
	

	# Organism information becomes note line
	# e.g. EntrezGeneId:Org:ID from upercase(Organism)/annotation
	if bcmlElement.find('Organism') is not None:
		hasNotes = True
		# Loop through all Organism entries
		for bcmlOrganism in bcmlElement.findall('Organism'):
			# Loop through all annotation entries
			for bcmlOrgAnnot in bcmlOrganism.findall('annotation'):
				# ENTREZ: ?
				organism = bcmlOrganism.attrib.get('name')
				if organism is None:
					organism = ""
				sbmlNotes += bcmlOrgAnnot.attrib.get('DB')+":"+organism.upper()+":"+bcmlOrgAnnot.attrib.get('ID').strip()+"\n"
	
	# Notes must be XHTML. !! For later: Store in Protein/gene notes
	
	# Return notes as one string
	return sbmlNotes

def idfy(string):
	str1 = re.sub("[-+():, ]", '_', string)
	if re.match('^[0-9]', str1):
		str1 = "_"+str1
	return str1


def main(argv):

	# Create an empty SBMLDocument object.  It's a good idea to check for
	# possible errors.  Even when the parameter values are hardwired like
	# this, it is still possible for a failure to occur (e.g., if the
	# operating system runs out of memory).
	try:
		document = SBMLDocument(2, 4)
	except ValueError:
		raise SystemExit('Could not create SBMLDocumention object')
	
	## SBML model
	
	# Create the basic Model object inside the SBMLDocument object.
	sbmlModel = document.createModel()
	# Check model correctly created
	check(sbmlModel, "create model")
	# Add a name to the model
	check(sbmlModel.setName(argv[1]), "Give name to model")
	
	# Set default units (best practice to set them)
	#check(sbmlModel.setTimeUnits("second"), 'set model-wide time units')
	#check(sbmlModel.setExtentUnits("mole"), 'set model units of extent')
	#check(sbmlModel.setSubstanceUnits('mole'), 'set model substance units')
	
	# Create species types which will be used for each Species creation
	#sbmlSTgene = SpeciesType(2, 4)
	#sbmlSTgene.setId("gene")
	#sbmlModel.addSpeciesType(sbmlSTgene)
	#sbmlSTrna = SpeciesType(2, 4)
	#sbmlSTrna.setId("RNA")
	#sbmlModel.addSpeciesType(sbmlSTrna)
	#sbmlSTprot = SpeciesType(2, 4)
	#sbmlSTprot.setId("protein")
	#sbmlModel.addSpeciesType(sbmlSTprot)
	#sbmlSTcomplex = SpeciesType(2, 4)
	#sbmlSTcomplex.setId("simple")
	#sbmlModel.addSpeciesType(sbmlSTcomplex)
	#sbmlSTcomplex = SpeciesType(2, 4)
	#sbmlSTcomplex.setId("complex")
	#sbmlModel.addSpeciesType(sbmlSTcomplex)
	
	## Read BCML and create elements in SBML model
	
	# Open BCML file and parse
	bcmlRoot = ET.parse(argv[1]).getroot()
	
	# Counter for compartments
	compartmentCounter = 1
	
	# Containers for Source, Sink, AndNote and OrNote
	sourceList = []
	sinkList = []
	andDict = {}
	orDict = {}
	
	# Loop through each compartment
	#print("* Compartment")
	for bcmlComp in bcmlRoot.iter('Compartment'):
		
		# Get Compartment name
		bcmlCompLabel = bcmlComp.attrib.get('label')
		bcmlCompLabel = re.sub('[\s+]', '', bcmlCompLabel)
		if bcmlCompLabel is None:
			bcmlCompLabel = 'default'
		#print(bcmlCompLabel)
		
		# Create an equivalent in SBML
		sbmlComp = sbmlModel.createCompartment()
		
		# Set id
		sbmlCompartmentId = "c"+str(compartmentCounter)
		check(sbmlComp.setId(sbmlCompartmentId), "set compartment Id")
		# Set name
		check(sbmlComp.setName(str(bcmlCompLabel)), "set compartment Name")
		# Set spatial dimensions (needed)
		check(sbmlComp.setSpatialDimensions(3), 'set compartment dimensions')
		# Set size (needed)
		check(sbmlComp.setSize(1), 'set compartment "size"')
		# Set units (needed)
		check(sbmlComp.setUnits("volume"), 'set compartment units')
		# Set outside
		#check(sbmlComp.setOutside("default"), 'set compartment outside')
		# Set constant (needed)
		check(sbmlComp.setConstant(True), 'set compartment constant')
				
		# Add macromolecules as species and related information as notes
		#print("* MacroMolecule")
		for bcmlMacromol in bcmlComp.findall('Macromolecule'):
			addMacroMolecule(bcmlMacromol, sbmlModel, sbmlCompartmentId)
		
		# Species:  NucleicAcidFeature : RNA, gene 
		# <UnitOfInformation label="gene"/"mRNA"
		#print("* Species")
		for bcmlNAfeature in bcmlComp.findall('NucleicAcidFeature'):
			addNucleicAcidFeature(bcmlNAfeature, sbmlModel, sbmlCompartmentId)
		
		# Species : SimpleChemical
		#print("* SimpleChemical")
		for bcmlSimpleChem in bcmlComp.findall('SimpleChemical'):
			addSimpleChemical(bcmlSimpleChem, sbmlModel, sbmlCompartmentId)
		
		# Species: Complex
		#print("* Complex")
		for bcmlComplex in bcmlComp.findall('Complex'):
			addComplex(bcmlComplex, sbmlModel, sbmlCompartmentId)
		
		# Species: Source, Sink
		# Are not explicit in SBML? eg reaction without reactant/product?
		#print("* Source / Sink")
		for bcmlSource in bcmlComp.findall('Source'):
			sourceList.append(bcmlSource.attrib.get('ID'))
		for bcmlSink in bcmlComp.findall('Sink'):
			sinkList.append(bcmlSink.attrib.get('ID'))

		# AndNode / OrNode
		#print("* AndNode / OrNode")
		for bcmlAndNode in bcmlComp.findall('AndNode'):
			andDict[bcmlAndNode.attrib.get('ID')] = [ idfy(str(bcmlLog.attrib.get('refNode'))) for bcmlLog in bcmlAndNode.findall('Logic') ]
		for bcmlOrNode in bcmlComp.findall('OrNode'):
			orDict[bcmlOrNode.attrib.get('ID')]	= [ idfy(str(bcmlLog.attrib.get('refNode'))) for bcmlLog in bcmlOrNode.findall('Logic') ]
		
		compartmentCounter+=1
	
	#print(andDict)
	#print(orDict)
		
	# Loop through each compartment again to create reactions
	# Needed because reactions add species references. If the corresponding species don't 
	# exist already, libsbml doesn't add them as reactant/product
	#print("* Compartment for reactions")
	sbmlReactionNb = 1
	for bcmlComp in bcmlRoot.iter('Compartment'):
		
		# Get Compartment name
		bcmlCompLabel = bcmlComp.attrib.get('label')
		bcmlCompLabel = re.sub('[\s+]', '', bcmlCompLabel)
		if bcmlCompLabel is None:
			bcmlCompLabel = 'default'
		#print(bcmlCompLabel)
		
		# Reaction
		# Process/Association/Dissociation
		#print("* Association")
		for bcmlReaction in bcmlComp.findall('Association'):
			addReaction(bcmlReaction, sbmlModel, sbmlCompartmentId, sbmlReactionNb, orDict, bcmlRoot)
			sbmlReactionNb += 1
		
		#print("* Dissociation")
		for bcmlReaction in bcmlComp.findall('Dissociation'):
			addReaction(bcmlReaction, sbmlModel, sbmlCompartmentId, sbmlReactionNb, orDict, bcmlRoot)
			sbmlReactionNb += 1
		
		#print("* Process")
		for bcmlReaction in bcmlComp.findall('Process'):
			addProcess(bcmlReaction, sbmlModel, sbmlCompartmentId, sbmlReactionNb, andDict, orDict, bcmlRoot)
			sbmlReactionNb += 1
	
	
	
	
	## Print SBML model in file
	
	# Print SBML in file
	outputdir = os.path.join(os.path.dirname(argv[1]), "to_SBML")
	if not os.path.exists(outputdir):
		os.makedirs(outputdir)
	outputfile = os.path.join(outputdir, os.path.splitext(os.path.basename(argv[1]))[0]+"_sbml.xml")
	writeSBMLToFile(document, outputfile)
	
	# Print SBML on STDOUT
	#print(writeSBMLToString(document))
	


if __name__ == "__main__":
	main(sys.argv)
