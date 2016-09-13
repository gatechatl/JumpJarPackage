package METABOLOMIC_DATABASE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.log4j.varia.ExternallyRolledFileAppender;

import CDKFunction.CalculateMonoisotope;
import CDKFunction.ConvertSMILE2Formula;
import METABOLOMICS.AIM.QueryMassWithPPM;
import MISC.ToolBox;

public class CreatePubchemDecoyDatabaseH {

	public static void execute(String[] args) {
		String fileName = args[0];			
		String outputFile = args[1];
		
		try {
			
			HashMap[] map = new HashMap[150000];
			for (int i = 0; i < 150000; i++) {
				map[i] = new HashMap();
			}
			HashMap all_formulas = new HashMap();
			
			int count = 0;
			
			FileWriter fwriter = new FileWriter(outputFile);
            BufferedWriter out = new BufferedWriter(fwriter);
			
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length == 7) {
					String smile = split[4];
					double raw_mass = Double.NaN; //new Double(split[6]);
					//if (ToolBox.check_formula_valid_element(split[1])) {
						if (!all_formulas.containsKey(split[1])) {
							double calc_mass = ToolBox.getMonoisotopicMass(split[1]);
							all_formulas.put(split[1], calc_mass);
							raw_mass = calc_mass;
						} else {
							raw_mass = (Double)all_formulas.get(split[1]);
						}
						int mass = (new Double(raw_mass * 100)).intValue();
						if (mass >= 150000) {
							mass = 149999;
						}
						count++;
						
						if (!map[mass].containsKey(split[1])) {
							
							String hillSystem = ToolBox.HillSystemOrder(split[1]);
							
							//String smileFormula = hillSystem;
							//boolean valid = false;
							//try {
							//	valid = QueryMassWithPPM.passRDBE(hillSystem);
							//} catch (Exception e) {
								//System.out.println(hillSystem + ": didn't work for RDBE");
							//}
							//if (valid) {
							//	if (raw_mass < 800) {
							//		if (smile.equals("")) {
							//			smileFormula = hillSystem;
							//		} else {
							//			try {
							//				smileFormula = ToolBox.HillSystemOrder(ConvertSMILE2Formula.SMILE2Formula(smile));
							//			} catch (Exception e) {
							//				smileFormula = hillSystem;
							//			}
							//		}
							//	}
							//	if (smileFormula.equals("")) {
							//		smileFormula = hillSystem;
							//	}
							//	//System.out.println(smileFormula);
							//	if (smileFormula.equals(hillSystem)) {
							//		System.out.println(count + "\t" + smileFormula + "\t" + smile);
									String hillSystem_DECOY = ToolBox.HillSystemOrder_DECOY_H(split[1]);					
									double decoy_mass = raw_mass + 1.0078250321;
									
									String text = hillSystem + ":" + split[6] + "|" + hillSystem_DECOY + ":" + decoy_mass;
									//map[mass].put(split[1], split[1]);
									map[mass].put(split[1], text);
									
									if (count % 10000 == 0) {
										//System.out.println(count);
										//out.flush();
									}
							//	} // end if smileFormula equals annotated formula
							//}
						}
					//}
				}
					
				
			}
			in.close();
			
			int max_length = 0;
			for (int i = 0; i < 150000; i++) {
				boolean first = true;
				Iterator itr = map[i].keySet().iterator();
				while (itr.hasNext()) {
					String formula = (String)itr.next();
					String str = (String)map[i].get(formula);
					if (first) {
						out.write(str);
					} else {
						out.write("\t" + str);
					}
					first = false;
				}
				out.write("\n");
				out.flush();
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
