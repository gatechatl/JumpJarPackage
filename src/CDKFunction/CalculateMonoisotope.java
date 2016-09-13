package CDKFunction;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import MISC.ToolBox;

public class CalculateMonoisotope {

	public static void main(String[] args) {
		double mass = calculateMonoisotopicPeakMass("C4H12N5Cl1");
		//double mass = calculateMonoisotopicPeakMass("C242H410N2O57");
		//mass = getMonoisotopicMassCDK("C");
		System.out.println(mass);
	}

	
	public static double calculateMonoisotopicPeakMass(String formula) {
		/*int C = ToolBox.retrieve_num_element(formula,  "C");
		int H = ToolBox.retrieve_num_element(formula,  "H");
		int N = ToolBox.retrieve_num_element(formula,  "N");
		int O = ToolBox.retrieve_num_element(formula,  "O");
		int P = ToolBox.retrieve_num_element(formula,  "P");
		int S = ToolBox.retrieve_num_element(formula,  "S");
		int Cl = ToolBox.retrieve_num_element(formula,  "Cl");
		int Br = ToolBox.retrieve_num_element(formula,  "Br");
		
		double mass = C * 12 + H * 1.0078250321 + N * 14.0030740052 + O * 15.9949146221 + P * 30.97376151 + S * 31.97207069;
		//return mass;*/
		return ToolBox.getMonoisotopicMass(formula);
	}
	public static IMolecularFormula getFormula(String str) {
		
		return MolecularFormulaManipulator.getMolecularFormula(str, DefaultChemObjectBuilder.getInstance());
	}
	public static IAtomContainer getAtom(String str) {
		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(str, DefaultChemObjectBuilder.getInstance());
		
		return MolecularFormulaManipulator.getAtomContainer(formula);
	}
	public static double getMonoisotopicMassBad(String str) {
		
		IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(str, DefaultChemObjectBuilder.getInstance());
		IAtomContainer atom = MolecularFormulaManipulator.getAtomContainer(formula);
		
		return AtomContainerManipulator.getNaturalExactMass(atom); 
	}
	/**
	 * This function can only calculate mass for molecules
	 *  
	 * @param str
	 * @return
	 */
	public static double getMonoisotopicMassCDK(String str) {
		return calculateMonoisotopicPeakMass(str);
	}
	
}
