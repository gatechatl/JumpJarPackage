package METABOLOMIC_DATABASE.PUBCHEM;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.Iterator;

import METABOLOMIC_DATABASE.PUBCHEM.DownloadPubchemXML.PUBCHEM;
import MISC.ToolBox;

/**
 * First generate a table for the processed PubchemFile
 * @author tshaw
 *
 */
public class ProcessPubChemXML {

	public static String parameter_info() {
		return "[xmlFile]";
	}
	public static void execute(String[] args) {
		
		try {
			
			//String inputFile = args[0];
			String inputPath = args[0];
			File f = new File(inputPath);
			for (File file: f.listFiles()) {
				
				String xmlFile = file.getPath();
				System.out.println("Reading file: " + xmlFile);
				if (!xmlFile.contains("gz")) {
					String outputFile = file.getPath() + ".txt";
					
					FileWriter fwriter = new FileWriter(outputFile);
		            BufferedWriter out = new BufferedWriter(fwriter);
					
					HashMap map = DownloadPubchemXML.readXML(xmlFile);
					Iterator itr = map.keySet().iterator();
					while (itr.hasNext()) {
						int id = (Integer)itr.next();
						PUBCHEM pubchem = (PUBCHEM)map.get(id);
						boolean radical = pubchem.RADICAL;
						boolean poscharge = pubchem.FORMULA.contains("+");
						boolean negcharge = pubchem.FORMULA.contains("-");
						boolean check_formula_valid = ToolBox.check_formula_valid_element(pubchem.FORMULA);
						boolean check_hydrogen_rule = ToolBox.check_hydrogen_rule(pubchem.FORMULA);
						String flag = radical + "\t" + poscharge + "\t" + negcharge + "\t" + check_formula_valid + "\t" + check_hydrogen_rule;
						out.write(pubchem.ID + "\t" + pubchem.FORMULA + "\t" + pubchem.InchIKey + "\t" + pubchem.InchI + "\t" + pubchem.SMILE + "\t" + pubchem.IUPAC_trad + "\t" + pubchem.MONOISOTOPICMASS + "\t" + flag + "\n");
					}
					out.close();
					file.delete();
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
