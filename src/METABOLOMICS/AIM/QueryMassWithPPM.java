package METABOLOMICS.AIM;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.List;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.interfaces.IMolecularFormula;
import org.openscience.cdk.tools.manipulator.MolecularFormulaManipulator;

import Interface.JumpInterface;
import MISC.ToolBox;
import MISC.RDBERule;
/**
 * Query a specific mass under a specified PPM
 * @author tshaw
 */
public class QueryMassWithPPM implements JumpInterface {
	public void execute(String[] args) {
		try {
			
			double queryMass = new Double(args[0]);
			double tolerance = new Double(args[1]);
			String folder = args[2];
			String rgdb = args[3];
			String term = "";
			if (args.length > 4) {
				term = args[4];
			}
			double diff = (tolerance * queryMass / 1E6);
			double lowerRange = queryMass - diff;
			double higherRange = queryMass + diff;
			
			int lower_index = (int)(lowerRange * 100);
			int higher_index = (int)(higherRange * 100);
			HashMap map = new HashMap();
			//LinkedList list = new LinkedList();
			for (int i = lower_index; i <= higher_index; i++) {
				File f = new File(folder + "/Formula" + (i + 1) + ".txt");
				if (f.isFile()) {
					String inputFile = folder + "/Formula" + (i + 1) + ".txt";
					FileInputStream fstream = new FileInputStream(inputFile);
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						String tag = "";
						HashMap tags = new HashMap();
						if (str.split("\t").length > 1) {
							tag = str.split("\t")[1];
							for (String tagstr: tag.split(",")) {
								tags.put(tagstr.toUpperCase(), tagstr);
							}
						}
						boolean search = false;
						
						if (tags.containsKey(term.toUpperCase()) || term.equals("")) {
							search = true;
						}
						if (search) {
							String[] split = str.split("\t")[0].split(":");						
							double mass = new Double(split[1]);
							if (mass >= lowerRange && mass <= higherRange) {
								String formula = split[0];
								String new_formula = ToolBox.standardize_name(formula);
								String new_str = new_formula + ":" + mass;
								if (rgdb.equals("yes")) {
									if (!(new_formula.equals("") || new_formula == null)) {
										//System.out.println("Check RDBE " + new_formula);
										if (passRDBE(new_formula)) {
											map.put(new_formula, "" + mass);
										}
									}
								} else if (rgdb.equals("no")) {
									if (!(new_formula.equals("") || new_formula == null)) {
										//System.out.println("Check RDBE " + new_formula);
										if (!passRDBE(new_formula)) {
											map.put(new_formula, "" + mass);
										}
									}
									
								} else {
									map.put(new_formula, "" + mass);
								}
								//if (!list.contains(new_str)) {
								//	list.add(new_str);
								//}
							} // if mass
						} // if serach
					}
					in.close();					
				}
			}
			Iterator itr = map.keySet().iterator();
			while (itr.hasNext()) {
				String key = (String)itr.next();
				String mass = (String)map.get(key);
				String outputKey = ToolBox.HillSystemOrder(key);
				System.out.println(outputKey + ":" + mass);
				//System.out.println(key + ":" + mass);
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	public static boolean passRDBE(String formulaStr) {
		RDBERule rdbeRule = new RDBERule();
		IMolecularFormula formula = getFormula(formulaStr);
		List<Double> list = rdbeRule.getRDBEValue(formula);
		if (list.size() > 0) {
			if (list.get(0) >= -1) {
				return true;
			}	
		}
		return false;		
	}
	public static IMolecularFormula getFormula(String str) {
		
		return MolecularFormulaManipulator.getMolecularFormula(str, DefaultChemObjectBuilder.getInstance());
	}
}
