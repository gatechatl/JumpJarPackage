package CDKFunction;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.smiles.SmilesGenerator;
import org.openscience.cdk.smiles.SmilesParser;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

public class ConvertSMILE2Formula {

	public static void main(String[] args) {
		
		try {
			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			//parse smiles
			IAtomContainer molecule = sp.parseSmiles("B(C(CC(C)C)NC(=O)C(CC1=CC=CC=C1)NC(=O)C2=NC=CN=C2)(O)O");
			IChemObjectBuilder builder = DefaultChemObjectBuilder.getInstance();;
			CDKHydrogenAdder.getInstance(builder).addImplicitHydrogens(molecule);
			SmilesGenerator sg = new SmilesGenerator();

			IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);
			
			System.out.println(MolecularFormulaManipulator.getString(formula));
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static String SMILE2Formula(String smile) {
		String result = "";	
		try {
			SmilesParser sp = new SmilesParser(DefaultChemObjectBuilder.getInstance());
			IAtomContainer molecule = sp.parseSmiles(smile);
			IMolecularFormula formula = MolecularFormulaManipulator.getMolecularFormula(molecule);
			
			result = MolecularFormulaManipulator.getString(formula);
		} catch (Exception e) {
			e.printStackTrace();
		}
		return result;
	}
	
}
