package METABOLOMICS.AIM;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.text.Collator;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.AtomContainerManipulator;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import CDKFunction.CalculateMonoisotope;
import Interface.JumpInterface;
import MISC.ToolBox;

/**
 * Given a formula, output all the available structures
 * @author tshaw
 *
 */
public class QueryStructureDatabase implements JumpInterface {

	public static void main(String[] args) {
		String formula = "C20H38FNO3Cl";
		double monoisotopicMass = CalculateMonoisotope.getMonoisotopicMassCDK(formula);
		System.out.println(monoisotopicMass);
	}
	
	
	public static String elementOnly(String formula) {
		LinkedList list = new LinkedList();
		String element = "";
		for (int i = 0; i < formula.length(); i++) {
			if (!isNumeric(formula.substring(i, i + 1))) {
				if (!element.equals("")) {
					if (Character.isUpperCase(formula.charAt(i))) {
						//System.out.println("Print something");
						//System.out.println(formula.charAt(i));
						list.add(element);
						element = formula.substring(i, i + 1);
					} else {
						element += formula.substring(i, i + 1);
					}
				} else {
					element += formula.substring(i, i + 1);
				}
				
			} else {
				
				list.add(element);
				element = "";
				
			}
		}
		if (!element.equals("")) {
			list.add(element);
		}
		sort(list);
		
		String result = "";
		Iterator itr = list.iterator();
		while (itr.hasNext()) {
			String str = (String)itr.next();
			result += str;
		}
		
		return result;
	}
	
	public static LinkedList sort(LinkedList list) {
		Collections.sort(list, new Comparator<String>() {
	         @Override
	         public int compare(String o1, String o2) {
	             return Collator.getInstance().compare(o1, o2);
	         }
	     });
		return list;
	}
	public static boolean isNumeric(String str) {
		return str.matches("-?\\d+(\\.\\d+)?");  //match a number with optional '-' and decimal.
	}
	
	@Override
	public void execute(String[] args) {
		try {
			
			String formula = args[0];
			String database = args[1];
			String datatype = args[2];
			String dbname = args[3];
			
			String elementOnly = elementOnly(formula);
			double monoisotopicMass = CalculateMonoisotope.getMonoisotopicMassCDK(formula);
			int mass = new Double(monoisotopicMass).intValue();
			//System.out.println(mass + "\t" + elementOnly);
			
			HashMap list = new HashMap();
			
			String fileName = database + "/" + mass + "/" + elementOnly + ".txt";
			File f = new File(database + "/" + mass + "/" + elementOnly + ".txt");
			//System.out.println(fileName);
			//System.out.println("Mass is this: " + mass);
			if (f.isFile()) {
				FileInputStream fstream = new FileInputStream(fileName);
				DataInputStream din = new DataInputStream(fstream);
				BufferedReader in = new BufferedReader(new InputStreamReader(din));
				while (in.ready()) {
					String str = in.readLine();
					String[] split = str.split("\t");
					
					if (split.length > 4) {
						String id = split[0];
						String database_formula = split[1];
						String inchikey = split[2];						
						String inchi = split[3];	
						String smile = "";
						String general_name = "";
						String dbname_tag = "";
						smile = split[4];
						
						String iupac = "";
						if (split.length > 5) {
							iupac = split[5];
						}
						//System.out.println(str);
						if (split.length >= 9) {
							general_name = split[7];
							dbname_tag = split[8];
							//System.out.println(general_name + "\t" + dbname_tag);
						}
						if (dbname_tag.equals("")) {
							dbname_tag = "pubchem";
						}
						//System.out.println(ToolBox.standardize_name(formula) + "\t" + ToolBox.standardize_name(database_formula));
						
						boolean grabrow = false;
						String[] dbname_split = dbname.split(",");
						for (String db: dbname_split) {
							if (dbname_tag.contains(db)) {
								grabrow = true;							
							}
						}
						/*if (dbname_tag.contains(dbname)) {
							grabrow = true;							
						}*/
						
						
						/*if (dbname.contains("PUBCHEM")) {
							grabrow = true;
							if (dbname_tag.contains("MISSINGPUBCHEM") || !inchi.contains("InChI=")) {
								grabrow = false;
							}
						}*/
						if (grabrow && ToolBox.standardize_name(formula).equals(ToolBox.standardize_name(database_formula))) {
							String output = database_formula;
							if (datatype.contains("ID")) {
								if (id.equals("")) {
									id = "NA";
								}
								output += "\t" + id;
							}
							if (datatype.contains("INCHIKEY")) {
								output += "\t" + inchikey;
							}
							if (datatype.contains("INCHISTR")) {
								output += "\t" + inchi;
							}
							if (datatype.contains("SMILE")) {
								output += "\t" + smile;
							}
							if (datatype.contains("IUPAC")) {
								if (iupac.equals("")) {
									iupac = "NA";
								}
								output += "\t" + iupac;
							}
							if (datatype.contains("GENERALNAME")) {
								if (general_name.equals("")) {
									general_name = "NA";
								}
								output += "\t" + general_name;
							}
							
							
							if (!list.containsKey(output)) {
								list.put(output, "");
								System.out.println(output);
							}
						}
					}
				}
				in.close();
			}
			/*File files = new File(database + "/" + monoisotopicMass);
			for (File file: files.listFiles()) {
			}*/
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
