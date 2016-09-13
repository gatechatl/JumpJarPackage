package MassDatabaseGeneration;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

import MISC.ToolBox;

public class CreateIndexFormulaExternalDB_H {
	public static void execute(String[] args) {
		try {
			HashMap all_formula = new HashMap();
			HashMap[] map = new HashMap[151000];
			for (int i = 0; i < 151000; i++) {
				map[i] = new HashMap();
			}

			String inputFile = args[0];
			String folder = args[1];
			double min = new Double(args[2]);
            double max = new Double(args[3]);

            String otherReference = args[4];
            String referenceNames = args[5];
            String formulaColIndexStr = args[6];
            String decoy_flags_str = args[7];
            String strict_valence_rule = args[8];
            boolean valence_rule_flag = false;
            if (strict_valence_rule.toUpperCase().equals("TRUE")) {
            	valence_rule_flag = true;
            } else {
            	valence_rule_flag = false;
            }
            int hashdivide = 20;
            HashMap[] sampleMap = new HashMap[hashdivide];
            for (int i = 0; i < hashdivide; i++) {
            	sampleMap[i] = new HashMap();
            }
            String[] formulaColIndexSplit = formulaColIndexStr.split(",");
            String[] splitRef = otherReference.split(",");
            String[] splitRefName = referenceNames.split(",");
            String[] splitDecoyFlag = decoy_flags_str.split(",");
            //HashMap[] maps = new HashMap[splitRef.length];
            for (int i = 0; i < splitRef.length; i++) {
            	int formulaColIndex = new Integer(formulaColIndexSplit[i]);
    			FileInputStream fstream = new FileInputStream(splitRef[i]);
    			DataInputStream din = new DataInputStream(fstream);
    			BufferedReader in = new BufferedReader(new InputStreamReader(din));
    			while (in.ready()) {
    				String str = in.readLine();
    				String[] split = str.split("\t");
    				
    				if (split.length > formulaColIndex) {
	    				String formula = split[formulaColIndex].split(":")[0].trim();
	    				if (ToolBox.check_formula_valid_element(formula)) {
	    					String name = "";
	    					double mass = ToolBox.getMonoisotopicMass(formula);
	    					if (splitDecoyFlag[i].toUpperCase().equals("D")) {
	    						//name = ToolBox.HillSystemOrder_DECOY_NH2(formula);
	    						name = ToolBox.HillSystemOrder_DECOY_H(formula);
	    						mass = ToolBox.getMonoisotopicMass(name);
	    						if (valence_rule_flag && ToolBox.check_hydrogen_rule(name)) {
	    							mass = -100;
	    						}
	    					} else {
	    						name = ToolBox.HillSystemOrder(formula);
	    						if (valence_rule_flag && !ToolBox.check_hydrogen_rule(name)) {
	    							mass = -100;
	    						}
	    					}
		    				
		    				if (min < mass && mass <= max) {
		    					double total = max - min;
		    					double part = total / hashdivide;
		    					for (int j = 0; j < hashdivide; j++) {
		    						double lower = min + part * j;
		    						double upper = min + part * (j + 1);
		    						if (lower < mass && mass <= upper) {
		    							if (sampleMap[j].containsKey(name)) {
		    								String stuff = (String)sampleMap[j].get(name);
		    								String[] listStuff = stuff.split(",");
		    								boolean foundExistDB = false;
		    								for (String db: listStuff) {
		    									if (db.equals(splitRefName[i])) {
		    										foundExistDB = true;
		    									}
		    								}
		    								if (!foundExistDB) {
		    									sampleMap[j].put(name, stuff + "," + splitRefName[i]);
		    								}
		    							} else {
		    								sampleMap[j].put(name, splitRefName[i]);
		    							}
		    						}    							
		    					}    					    					
		    				}
	    				}
    				}
    			}
    			in.close();    			
            }
            
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length > 2) {
					for (int i = 2; i < split.length; i++) {
						String[] split2 = split[i].split(":");
						String formula = split2[0];
						double mass = new Double(split2[1]);					
						int index = (int)(mass * 100);
						if (mass >= min && mass < max) {
							formula = ToolBox.HillSystemOrder(formula);
							map[index].put(formula + ":" + mass, "");
						}
					}
				}
			}
			in.close();			
			
			File f = new File(folder);
			if (!f.exists()) {
				f.mkdir();
			}
			
			for (int i = 0; i < 151000; i++) {
				if (i >= min * 100 && i < max * 100) {
					FileWriter fwriter = new FileWriter(folder + "/Formula" + (i + 1) + ".txt"); 
					BufferedWriter out = new BufferedWriter(fwriter);
					
					Iterator itr = map[i].keySet().iterator();
					while (itr.hasNext()) {
						String str = (String)itr.next();
						String[] split = str.split(":");
						String formula = ToolBox.HillSystemOrder(split[0]);
						boolean pass_valence_rule = false;
						if (ToolBox.check_hydrogen_rule(formula)) {
							pass_valence_rule = true;
						} else {
							pass_valence_rule = false;
						}
						double mass = ToolBox.getMonoisotopicMass(formula);
    					double total = max - min;
    					double part = total / hashdivide;
    					String databaseName = "";
    					for (int j = 0; j < hashdivide; j++) {
    						double lower = min + part * j;
    						double upper = min + part * (j + 1);
    						if (lower < mass && mass <= upper) {
    							if (sampleMap[j].containsKey(formula)) {
    								databaseName = (String)sampleMap[j].get(formula);
    							}
    						}
    					}
    					if (databaseName.equals("")) {
    						if (pass_valence_rule) {
    							databaseName = "V_target";
    						} else {
    							databaseName = "V_decoy";
    						}
    					} else {
    						if (pass_valence_rule) {
    							databaseName += ",V_target";
    						} else {
    							databaseName += ",V_decoy";
    						}
    					}
    						
    					
    					if (databaseName.equals("")) {
    						
    						out.write(str + "\n");
    					} else {
    						out.write(str + "\t" + databaseName + "\n");
    					}
						/*String[] split = str.split(":");
						String formula = ToolBox.standardize_name(split[0]);
						String mass = split[1];
						String hits = "";
						for (int j = 0; j < maps.length; j++) {
							if (maps[j].containsKey(formula)) {
								hits += splitRefName[j] + ",";
							}
						}
						if (hits.equals("")) {
							out.write(formula + ":" + mass + "\n");
						} else {
							out.write(formula + ":" + mass + "\t" + hits + "\n");
						}*/
					}
					out.close();
				}
			}						
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

