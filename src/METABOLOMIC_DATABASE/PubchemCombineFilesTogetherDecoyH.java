package METABOLOMIC_DATABASE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.HashMap;
import java.util.Iterator;

import CDKFunction.CalculateMonoisotope;
import CDKFunction.ConvertSMILE2Formula;
import METABOLOMICS.AIM.QueryMassWithPPM;
import MISC.ToolBox;

public class PubchemCombineFilesTogetherDecoyH {
	public static void execute(String[] args) {

		try {
			File dir = new File(args[0]);
			
			FileWriter fwriter = new FileWriter(args[1]); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Formula_Expansion");
			BufferedWriter out = new BufferedWriter(fwriter);
			
			HashMap all_formulas = new HashMap();
			HashMap[] map = new HashMap[150000];
			for (int i = 0; i < 150000; i++) {
				map[i] = new HashMap();
			}
			
			double count = 0;
			for (File f: dir.listFiles()) {
				if (f.isDirectory()) {
					for (File d: f.listFiles()) {
						if (!d.getName().contains("+") && !d.getName().contains("-")) {
							String inputFile = d.getPath();
							FileInputStream fstream = new FileInputStream(inputFile);
							DataInputStream din = new DataInputStream(fstream);
							BufferedReader in = new BufferedReader(new InputStreamReader(din));
							while (in.ready()) {
								String str = in.readLine();
								
								out.write(str + "\n");
								String[] split = str.split("\t");
								if (split.length == 7) {
									String smile = split[4];
									if (ToolBox.check_formula_valid_element(split[1])) {
										double raw_mass = Double.NaN; //new Double(split[6]);
										
										if (!all_formulas.containsKey(split[1])) {
											//double calc_mass = CalculateMonoisotope.getMonoisotopicMassCDK(split[1]);
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
											/*boolean valid = false;
											try {
												valid = QueryMassWithPPM.passRDBE(hillSystem);
											} catch (Exception e) {
												//System.out.println(hillSystem + ": didn't work for RDBE");
											}*/
											//if (valid) {
												/*if (raw_mass < 1500) {
													if (smile.equals("")) {
														smileFormula = hillSystem;
													} else {
														try {
															smileFormula = ToolBox.HillSystemOrder(ConvertSMILE2Formula.SMILE2Formula(smile));
														} catch (Exception e) {
															smileFormula = hillSystem;
														}
													}
												}
												if (smileFormula.equals("")) {
													smileFormula = hillSystem;
												}*/
												//System.out.println(smileFormula);
												//if (smileFormula.equals(hillSystem)) {
													//System.out.println(count + "\t" + smileFormula + "\t" + smile);
													String hillSystem_DECOY = ToolBox.HillSystemOrder_DECOY_H(split[1]);					
													double decoy_mass = raw_mass + 1.0078250321;
													
													double calc_decoy_mass = ToolBox.getMonoisotopicMass(hillSystem_DECOY);
													//String text = hillSystem + ":" + split[6] + "|" + hillSystem_DECOY + ":" + calc_decoy_mass;
													String text = hillSystem + ":" + raw_mass + "|" + hillSystem_DECOY + ":" + calc_decoy_mass;
													//map[mass].put(split[1], split[1]);
													map[mass].put(split[1], text);
													
													if (count % 10000 == 0) {
														//System.out.println(count);
														//out.flush();
													}
												//} // end if smileFormula equals annotated formula
											//}
										}
										/*if (!map[mass].containsKey(split[1])) {
											
											String hillSystem = ToolBox.HillSystemOrder(split[1]);
											String hillSystem_DECOY = ToolBox.HillSystemOrder_DECOY(split[1]);					
											double decoy_mass = raw_mass + 1.0078250321;
											
											String text = hillSystem + ":" + raw_mass + "|" + hillSystem_DECOY + ":" + decoy_mass;
											//map[mass].put(split[1], split[1]);
											map[mass].put(split[1], text);
											
											
										}*/
									} // end if toolbox
								} // end split length == 7
								count++;
							}
							in.close();
						}
					}
					out.flush();
				}
			}
			out.close();
			System.out.println(count);
			
			FileWriter fwriter2 = new FileWriter(args[2]); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Formula_Expansion");
			BufferedWriter out2 = new BufferedWriter(fwriter2);
			
			FileWriter fwriter3 = new FileWriter(args[3]); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\Formula_Expansion");
			BufferedWriter out3 = new BufferedWriter(fwriter3);
			
			int max_length = 0;
			for (int i = 0; i < 150000; i++) {
				boolean first = true;
				Iterator itr = map[i].keySet().iterator();
				while (itr.hasNext()) {
					String formula = (String)itr.next();
					String str = (String)map[i].get(formula);
					out3.write(str.replaceAll("\\|", "\t") + "\n");
					if (first) {
						out2.write(i + "\t" + str);
					} else {
						out2.write("\t" + str);
					}
					first = false;
				}
				if (!first) {
					out2.write("\n");
				}
				out2.flush();
			}
			out2.close();
			out3.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

