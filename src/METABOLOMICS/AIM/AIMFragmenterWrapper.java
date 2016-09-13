package METABOLOMICS.AIM;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.InputStreamReader;
import java.util.UUID;
import MISC.CommandLine;

public class AIMFragmenterWrapper {

	public static String parameter_info() {
		return "[jumpjarPath] [tmpPath] [smile string] [mass (cutoff)] [depth (number)] [aromatic_ring_flag yes/no] [neutralLossFile] [bondEnergyFile]";
	}
	public static void execute(String[] args) {
		
		try {
			String jumpjarPath = args[0];
			String tmpPath = args[1];
			String smile = args[2]; //"C1=CC(=C(C=C1C2=C(C(=O)C3=C(C=C(C=C3O2)O)O)OC4C(C(C(C(O4)CO)O)O)O)O)O";
			
			double mass = new Double(args[3]); //303.0499;
			
			int depth = new Integer(args[4]);
			
			boolean aromatic_ring_flag = false;
			String aromatic_ring_str = args[5];
			if (aromatic_ring_str.equals("yes")) {
				aromatic_ring_flag = true;
			}
			String molecularFormulaRedundancyCheck_str = args[6];
			boolean molecularFormulaRedundancyCheck = true;
			if (molecularFormulaRedundancyCheck_str.equals("no")) {
				molecularFormulaRedundancyCheck = false;
			}
			String neutralLossFile = args[7];	
			String bondEnergyFile = args[8];
			
			String buffer = UUID.randomUUID().toString();
			File f = new File(tmpPath);
			if (!f.exists()) {
				f.mkdir();
			}
			String tmpOutput = tmpPath + "/" + buffer;
			//"[jumpjarPath] [smile string] [mass (cutoff)] [depth (number)] [aromatic_ring_flag yes/no] [neutralLossFile] [bondEnergyFile]";
			String command = jumpjarPath + " -FragmentSMILE \"" + smile + "\" " + mass + " " + depth + " " + aromatic_ring_str + " " + molecularFormulaRedundancyCheck_str + " " + neutralLossFile + " " + bondEnergyFile + " > " +  tmpOutput;
			
			//System.out.println(command);
			CommandLine.executeCommand(command);
			
			//System.out.println("Reading tmpOutput:");
			f = new File(tmpOutput);
			if (f.exists()) {
				boolean error = false;
				FileInputStream fstream = new FileInputStream(tmpOutput);
				DataInputStream din = new DataInputStream(fstream);
				BufferedReader in = new BufferedReader(new InputStreamReader(din));
				while (in.ready()) {
					String str = in.readLine();
					if (str.contains("Error!")) {
						error = true;
					}										
				}
				in.close();
				
				if (error) {
					System.out.println("Found Error");
				} else {
					fstream = new FileInputStream(tmpOutput);
					din = new DataInputStream(fstream);
					in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						System.out.println(str);
					}
					in.close();
				}
			}			
			f.delete();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
