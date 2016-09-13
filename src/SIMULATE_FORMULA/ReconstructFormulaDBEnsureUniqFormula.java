package SIMULATE_FORMULA;

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

public class ReconstructFormulaDBEnsureUniqFormula {

	public static void execute(String[] args) {
		
		try {
			String filePath = args[0];
			String outputPath = args[1];
			
			File folder = new File(filePath);
			
			int count = 0;
			File[] listOfFiles = folder.listFiles();
			for (int i = 0; i < listOfFiles.length; i++) {
				count++;
				System.out.println(count);
				if (listOfFiles[i].isFile()) {
					
					HashMap map = new HashMap();
					FileInputStream fstream = new FileInputStream(listOfFiles[i].getPath());
					DataInputStream din = new DataInputStream(fstream);
					BufferedReader in = new BufferedReader(new InputStreamReader(din));
					while (in.ready()) {
						String str = in.readLine();
						String formula = ToolBox.HillSystemOrder(str.split(":")[0]);
						map.put(formula, str.split(":")[1]);
					}
					in.close();
					
					FileWriter fwriter = new FileWriter(outputPath + "/" + listOfFiles[i].getName()); //"C:\\Users\\tshaw\\Desktop\\METABOLOMICS\\MISSILE\\All_Formulas.txt");
					BufferedWriter out = new BufferedWriter(fwriter);
					Iterator itr = map.keySet().iterator();
					while (itr.hasNext()) {
						String formula = (String)itr.next();
						String mass = (String)map.get(formula);
						out.write(formula + ":" + mass + "\n");
					}
					out.close();
				}
			}
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
