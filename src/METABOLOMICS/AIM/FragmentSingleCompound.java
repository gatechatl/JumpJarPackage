package METABOLOMICS.AIM;

import java.io.File;
import java.util.ArrayList;
import java.util.List;
import java.util.Vector;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.interfaces.IMolecule;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;

import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import fragmenter.Fragmenter;
import de.ipbhalle.metfrag.massbankParser.Peak;
import de.ipbhalle.metfrag.read.Molfile;
import de.ipbhalle.metfrag.spectrum.WrapperSpectrum;
import de.ipbhalle.metfrag.tools.MolecularFormulaTools;
import de.ipbhalle.metfrag.tools.Render;

public class FragmentSingleCompound {
	
	private int treeDepth;
	private boolean sumFormulaRedundancyCheck;
	private List<Double> dissociationEnergyToReturn;
	private List<String> fromNeutralLoss;
	
	/**
	 * Instantiates a new fragment single compound with default options:
	 * <ul>
	 * 	<li> tree depth = 2
	 * 	<li> sum formula redundancy check = true
	 *  <li> brak aromatic rings = true
	 * </ul>
	 */
	public FragmentSingleCompound()
	{
		setTreeDepth(2);
		setSumFormulaRedundancyCheck(true);
		dissociationEnergyToReturn = new ArrayList<Double>();
		fromNeutralLoss = new ArrayList<String>();
	}
	
	/**
	 * Gets the fragments.
	 * 
	 * @param smilesToFragment the smiles to fragment
	 * @param minMass the min mass
	 * @param treeDepth the tree depth
	 * @param Render the render
	 * 
	 * @return the fragments
	 * @throws Exception 
	 */
	public List<String> getFragments(String smilesToFragment, Double minMass, boolean render, boolean aromatic_ring_flag) throws Exception
	{
		List<String> resultingSmiles = new ArrayList<String>();
		SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
		//parse smiles
		IAtomContainer molecule = sp.parseSmiles(smilesToFragment);
		//configure atoms
		AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
		//add all hydrogens explicitly
		CDKHydrogenAdder adder1 = CDKHydrogenAdder.getInstance(molecule.getBuilder());
        adder1.addImplicitHydrogens(molecule);
        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule); 
        
        Double molMass = MolecularFormulaTools.getMonoisotopicMass(MolecularFormulaManipulator.getMolecularFormula(molecule));
		molMass = (double)Math.round((molMass)*10000)/10000;
        
		//this creates a pseudo peak list --> only fragments heavier than this peak will be generated
		//just a dirty hack...
		String peakString = minMass.toString() + " 10000 999\n"; 
		WrapperSpectrum spectrum = new WrapperSpectrum(peakString, 1, molMass);		
		
		//constructor for fragmenter
		//Fragmenter fragmenter = new Fragmenter((Vector<Peak>)spectrum.getPeakList().clone(), minMass,  true, sumFormulaRedundancyCheck, false);
		Fragmenter fragmenter = new Fragmenter((Vector<Peak>)spectrum.getPeakList().clone(), minMass,  aromatic_ring_flag, sumFormulaRedundancyCheck, false);
		List<IAtomContainer> listOfFrags = fragmenter.generateFragmentsInMemory(molecule, false, this.treeDepth);
		
		List<String> results = new ArrayList<String>();
		
		for (IAtomContainer fragment : listOfFrags) {
			SmilesGenerator sg = new SmilesGenerator();
			sg.setUseAromaticityFlag(true);
			IMolecule mol = new Molecule(fragment);
			results.add(sg.createSMILES(mol));
			if(fragment.getProperty("BondEnergy") == null)
				dissociationEnergyToReturn.add(0.0);
			else
				dissociationEnergyToReturn.add(Double.parseDouble((String)fragment.getProperty("BondEnergy")));
			
			if(fragment.getProperty("NeutralLossRule") == null)
				fromNeutralLoss.add("");
			else
				fromNeutralLoss.add((String)fragment.getProperty("NeutralLossRule"));
			
//			System.out.println(fragment.getProperty((String)fragment.getProperty("BondEnergy")));
		}
		
		if(render)
			Render.Draw(molecule, listOfFrags, "Fragments");
		
		return results;
	}
	
	/**
	 * Gets the energies in the same order as the fragments were returned.
	 * 
	 * @return the energies
	 */
	public List<Double> getEnergies()
	{
		return dissociationEnergyToReturn;
	}
	
	
	/**
	 * Gets the neutral losses.
	 * 
	 * @return the neutral losses
	 */
	public List<String> getNeutralLosses()
	{
		return fromNeutralLoss;
	}
	
	/**
	 * Sets the tree depth.
	 * 
	 * @param treeDepth the new tree depth
	 */
	public void setTreeDepth(int treeDepth) {
		this.treeDepth = treeDepth;
	}

	/**
	 * Gets the tree depth.
	 * 
	 * @return the tree depth
	 */
	public int getTreeDepth() {
		return treeDepth;
	}

	/**
	 * Sets the sum formula redundancy check.
	 * 
	 * @param sumFormulaRedundancyCheck the new sum formula redundancy check
	 */
	public void setSumFormulaRedundancyCheck(boolean sumFormulaRedundancyCheck) {
		this.sumFormulaRedundancyCheck = sumFormulaRedundancyCheck;
	}

	/**
	 * Checks if is sum formula redundancy check.
	 * 
	 * @return true, if is sum formula redundancy check
	 */
	public boolean isSumFormulaRedundancyCheck() {
		return sumFormulaRedundancyCheck;
	}
	
	public static void calculate(String smile, double mass, int depth, boolean aromatic_ring_flag, boolean molecularFormulaRedundancyCheck) {
		//example values
		//String smiles = "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)CO)O)O)O)O)O";
		String smiles = smile;
		//Double minMass = 303.0499;
		Double minMass = mass;
		Boolean render = false;
		//get command line arguments
		
		FragmentSingleCompound test = new FragmentSingleCompound();
		//those are the default values...no need to set them like this
		//test.setSumFormulaRedundancyCheck(true);
		test.setSumFormulaRedundancyCheck(molecularFormulaRedundancyCheck);
		test.setTreeDepth(depth);
		
		try {
			List<String> resultingFragments = test.getFragments(smiles, minMass, render, aromatic_ring_flag);
			List<Double> resultingEnergies = test.getEnergies();
			List<String> resultingNeutralLosses = test.getNeutralLosses();
			
			System.out.println("Fragment count: " + resultingFragments.size());
			for (int i = 0; i < resultingFragments.size(); i++) {
				
				SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
				//parse smiles
				IAtomContainer molecule = sp.parseSmiles(resultingFragments.get(i));
				//configure atoms
//				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//				//add all hydrogens explicitly
//				CDKHydrogenAdder adder1 = CDKHydrogenAdder.getInstance(molecule.getBuilder());
//		        adder1.addImplicitHydrogens(molecule);
//		        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule); 
		        
		        IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);		  			
				
				System.out.print(resultingFragments.get(i));
				System.out.print(" " + resultingEnergies.get(i) + " " + MolecularFormulaManipulator.getString(formula) + " " + MolecularFormulaTools.getMonoisotopicMass(formula) + " " + resultingNeutralLosses.get(i) +"\n");
			}
		} catch (Exception e) {
			System.out.println("Error! TODO...");
			e.printStackTrace();
		}
		
	}
	
	public static void main(String[] args) {
		//example values
		String smiles = "C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)CO)O)O)O)O)O";
		Double minMass = 303.0499;
		Boolean render = false;
		boolean aromatic_ring_flag = false;
		//get command line arguments
		if(args != null && args.length == 3)
		{
//			smiles = args[0];
//			minMass = Double.parseDouble(args[1]);
			if(args[2].equals("1"))
				render = true;
		}
		else
		{
			System.err.println("Please enter CL values!\n1. value: Smiles to fragment\n2. value: Minimum mass\n3. value: Render fragments? (1 --> true, 0 --> false)\nExample: C1C(OC2=CC(=CC(=C2C1=O)O)O)C3=CC=C(C=C3)O 42.1 1\n");
			System.exit(1);
		}
		
		FragmentSingleCompound test = new FragmentSingleCompound();
		//those are the default values...no need to set them like this
		test.setSumFormulaRedundancyCheck(true);
		test.setTreeDepth(2);
		
		try {
			List<String> resultingFragments = test.getFragments(smiles, minMass, render, aromatic_ring_flag);
			List<Double> resultingEnergies = test.getEnergies();
			List<String> resultingNeutralLosses = test.getNeutralLosses();
			
			System.out.println("Fragment count: " + resultingFragments.size());
			for (int i = 0; i < resultingFragments.size(); i++) {
				
				SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
				//parse smiles
				IAtomContainer molecule = sp.parseSmiles(resultingFragments.get(i));
				//configure atoms
//				AtomContainerManipulator.percieveAtomTypesAndConfigureAtoms(molecule);
//				//add all hydrogens explicitly
//				CDKHydrogenAdder adder1 = CDKHydrogenAdder.getInstance(molecule.getBuilder());
//		        adder1.addImplicitHydrogens(molecule);
//		        AtomContainerManipulator.convertImplicitToExplicitHydrogens(molecule); 
		        
		        IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);
		  			
				
				System.out.print(resultingFragments.get(i));
				System.out.print(" " + resultingEnergies.get(i) + " " + MolecularFormulaManipulator.getString(formula) + " " + MolecularFormulaTools.getMonoisotopicMass(formula) + " " + resultingNeutralLosses.get(i) +"\n");
			}
		} catch (Exception e) {
			System.out.println("Error! TODO...");
			e.printStackTrace();
		}
		
	}

}
