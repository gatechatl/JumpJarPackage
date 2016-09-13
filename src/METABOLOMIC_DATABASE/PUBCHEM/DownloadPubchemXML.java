package METABOLOMIC_DATABASE.PUBCHEM;


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.net.URL;
import java.nio.channels.Channels;
import java.nio.channels.ReadableByteChannel;
import java.util.HashMap;
import java.util.Iterator;
import java.util.UUID;


import MISC.ToolBox;

public class DownloadPubchemXML {

	public static void execute(String[] args) {
		
		try {
			String outputDirectory = args[0];			
			String outputFile = args[1];
			
			File file = new File(outputDirectory);
			if (!file.exists()) {
				file.mkdir();
			}
			FileWriter fwriter = new FileWriter(outputFile);
            BufferedWriter out = new BufferedWriter(fwriter);
			
            FileWriter fwriter_log = new FileWriter("DownloadPubchemXML.log");
            BufferedWriter out_log = new BufferedWriter(fwriter_log);
			
            double bad_entries = 0;
            double good_entries = 0;
            
			int start = 1;
			int end = 25000;
			while (end < 73525000) {
				String fileName = "Compound_" + convertInt2String(start) + "_" + convertInt2String(end) + ".xml.gz";
				
				
				String outputFile2 = outputDirectory + "/Compound_" + convertInt2String(start) + "_" + convertInt2String(end) + ".summary.txt";;
	        	
				out.write(outputFile2 + "\n");
				out.flush();
				FileWriter fwriter2 = new FileWriter(outputFile2);
	            BufferedWriter out2 = new BufferedWriter(fwriter2);
				
				String command = "wget ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/XML/" + fileName;
				
				//System.out.println(command);
				downloadFile("ftp://ftp.ncbi.nih.gov/pubchem/Compound/CURRENT-Full/XML/" + fileName, fileName);
				//executeCommand(command);
				start = start + 25000;
				end = end + 25000;				
				command = "gunzip " + fileName;
				executeCommand(command);
				
				

				String xmlFile = fileName.replaceAll("\\.gz", "");;
				HashMap map = readXML(xmlFile);
				Iterator itr = map.keySet().iterator();
				while (itr.hasNext()) {
					int id = (Integer)itr.next();
					PUBCHEM pubchem = (PUBCHEM)map.get(id);
					
					if (!pubchem.RADICAL && !pubchem.FORMULA.contains("+") && !pubchem.FORMULA.contains("-")) {
						if (ToolBox.check_formula_valid_element(pubchem.FORMULA) && ToolBox.check_hydrogen_rule(pubchem.FORMULA)) {
							out2.write(pubchem.ID + "\t" + pubchem.FORMULA + "\t" + pubchem.InchIKey + "\t" + pubchem.InchI + "\t" + pubchem.SMILE + "\t" + pubchem.IUPAC_trad + "\t" + pubchem.MONOISOTOPICMASS + "\n");
							out2.flush();
							good_entries++;
						} else {
							bad_entries++;
						}
					} else {
						bad_entries++;
					}
					//System.out.println(pubchem.ID + "\t" + pubchem.FORMULA + "\t" + pubchem.InchIKey + "\t" + pubchem.InchI + "\t" + pubchem.SMILE + "\t" + pubchem.IUPAC_trad + "\t" + pubchem.MONOISOTOPICMASS);
				}
				out2.close();
				
				command = "rm -rf " + xmlFile;
				executeCommand(command);
				
			}
			
			out.close();
			
			out_log.write("Good Entries: " + good_entries + "\n");
			out_log.write("Bad Entries: " + bad_entries + "\n");
			out_log.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}
	
	public static void downloadFile(String path, String fileName) {
		try {
			URL website = new URL(path);
			ReadableByteChannel rbc = Channels.newChannel(website.openStream());
			FileOutputStream fos = new FileOutputStream(fileName);
			fos.getChannel().transferFrom(rbc, 0, Long.MAX_VALUE);
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	public static HashMap readXML(String fileName) {
		HashMap map = new HashMap();
		try {
			PUBCHEM pubchem = new PUBCHEM();
			int id = 0;
			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				if (str.contains("PC-CompoundType_id_cid")) {
					if (!pubchem.ID.equals("")) {
						map.put(id, pubchem);
						//System.out.println("Reading ID: " + id);
					}
					String idstr = str.replaceAll("PC-CompoundType_id_cid", "");
					idstr = idstr.replaceAll("<>", "");
					idstr = idstr.replaceAll("</>", "").trim();
					id = new Integer(idstr);					
					pubchem = new PUBCHEM();
					pubchem.ID = idstr;
				}
				
				if (str.contains("Molecular Formula")) {
					while (!str.contains("PC-InfoData_value_sval")) {
						str = in.readLine();
					}
					if (str.contains("PC-InfoData_value_sval")) {
						String formula = str.replaceAll("PC-InfoData_value_sval", "");
						formula = formula.replaceAll("<>", "");
						formula = formula.replaceAll("</>", "").trim();
						//System.out.println(formula);
						pubchem.FORMULA = formula;
					}
				}
				if (str.contains("<PC-Atoms_radical>")) {
					pubchem.RADICAL = true;
				}
				if (str.contains("SMILES")) {
					str = in.readLine();
					if (str.contains("Canonical")) {
						while (!str.contains("PC-InfoData_value_sval")) {
							str = in.readLine();
						}
						if (str.contains("PC-InfoData_value_sval")) {
							String val = str.replaceAll("PC-InfoData_value_sval", "");
							val = val.replaceAll("<>", "");
							val = val.replaceAll("</>", "").trim();
							//System.out.println(val);
							pubchem.SMILE = val;
						}
					}
				}
				
				if (str.contains("<PC-Urn_label>InChI</PC-Urn_label>")) {
					while (!str.contains("PC-InfoData_value_sval")) {
						str = in.readLine();
					}
					if (str.contains("PC-InfoData_value_sval")) {
						String val = str.replaceAll("PC-InfoData_value_sval", "");
						val = val.replaceAll("<>", "");
						val = val.replaceAll("</>", "").trim();
						//System.out.println(val);
						pubchem.InchI = val;
					}
				}
				
				if (str.contains("<PC-Urn_label>InChIKey</PC-Urn_label>")) {
					while (!str.contains("PC-InfoData_value_sval")) {
						str = in.readLine();
					}
					if (str.contains("PC-InfoData_value_sval")) {
						String val = str.replaceAll("PC-InfoData_value_sval", "");
						val = val.replaceAll("<>", "");
						val = val.replaceAll("</>", "").trim();
						//System.out.println(val);
						pubchem.InchIKey = val;
					}
				}

				if (str.contains("<PC-Urn_name>MonoIsotopic</PC-Urn_name>")) {
					while (!str.contains("PC-InfoData_value_fval")) {
						str = in.readLine();
					}
					if (str.contains("PC-InfoData_value_fval")) {
						String val = str.replaceAll("PC-InfoData_value_fval", "");
						val = val.replaceAll("<>", "");
						val = val.replaceAll("</>", "").trim();
						//System.out.println(val);
						pubchem.MONOISOTOPICMASS = val;
					}
				}
				if (str.contains("<PC-Urn_label>Log P</PC-Urn_label>")) {
					while (!str.contains("PC-InfoData_value_fval")) {
						str = in.readLine();
					}
					if (str.contains("PC-InfoData_value_fval")) {
						String val = str.replaceAll("PC-InfoData_value_fval", "");
						val = val.replaceAll("<>", "");
						val = val.replaceAll("</>", "").trim();
						//System.out.println(val);
						pubchem.LogP = val;
					}
				}
				
				
				if (str.contains("IUPAC Name")) {
					str = in.readLine();
					if (str.contains("Traditional")) {
						while (!str.contains("PC-InfoData_value_sval")) {
							str = in.readLine();
						}
						if (str.contains("PC-InfoData_value_sval")) {
							String val = str.replaceAll("PC-InfoData_value_sval", "");
							val = val.replaceAll("<>", "");
							val = val.replaceAll("</>", "").trim();
							//System.out.println(val);
							pubchem.IUPAC_trad = val;
						}
					}
				}
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
		return map;
	}
	public static String convertInt2String(int num) {
		String str = (new Integer(num)).toString();
		while (str.length() < 9) {
			str = "0" + str;
		}
		return str;
	}
	public static void executeCommand(String executeThis) {
    	try {
    		
    		String buffer = UUID.randomUUID().toString();
        	writeFile(buffer + "tempexecuteCommand.sh", executeThis);
        	String[] command = {"sh", buffer + "tempexecuteCommand.sh"};
        	Process p1 = Runtime.getRuntime().exec(command);
        	BufferedReader inputn = new BufferedReader(new InputStreamReader(p1.getInputStream()));
        	String line=null;
        	while((line=inputn.readLine()) != null) {}
        	inputn.close();
        	p1.destroy();
        	File f = new File(buffer + "tempexecuteCommand.sh");
        	f.delete();
        } catch (Exception e) {
        	e.printStackTrace();
        }
	}
	public static void writeFile(String fileName, String command) {
    	try {
        	FileWriter fwriter2 = new FileWriter(fileName);
            BufferedWriter out2 = new BufferedWriter(fwriter2);
            out2.write(command + "\n");
            out2.close();
        } catch (Exception e) {
        	e.printStackTrace();
        }
	}

	static class PUBCHEM {
		public String ID = "";
		public String FORMULA = "";
		public String LogP = "";
		public String InchI = "";
		public String InchIKey = "";
		public String SMILE = "";
		public String IUPAC_trad= "";
		public String MONOISOTOPICMASS = "";
		public boolean RADICAL = false;	
	}
}

