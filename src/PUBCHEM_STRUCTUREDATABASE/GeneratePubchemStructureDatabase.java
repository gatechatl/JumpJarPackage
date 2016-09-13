package PUBCHEM_STRUCTUREDATABASE;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.text.Collator;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.LinkedList;
import java.util.UUID;

import MISC.ToolBox;

public class GeneratePubchemStructureDatabase {

	public static String parameter_info() {
		return "[fileName] [outputFolder]";
	}
	public static void execute(String[] args) {

		try {
			
			int count = 0;
			HashMap complete = new HashMap();
			//LinkedList fileList = file2List(args[0]);
			//Iterator itr = fileList.iterator();
			String fileName = args[0];
			String outputFolder = args[1];

			for (int i = 1; i < 10000; i++) {
				//File massDir = new File("Indexed/" + i);
				File massDir = new File(outputFolder + "/" + i);
				if (!massDir.exists()) {
					massDir.mkdir();
				}
			}
			

			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String id = split[0];
				String formula = split[1];
				String inchikey = split[2];
				String inchi = split[3];
				String smile = split[4];
				String iupac = split[5];
				
				
				//if (ToolBox.check_formula_valid_element(split[1])) {
					double calc_mass = ToolBox.getMonoisotopicMass(split[1]);
					String monomass = calc_mass + "";
					/*
					PUBCHEM pubchem = new PUBCHEM();
					pubchem.ID = id;
					pubchem.FORMULA = formula;
					pubchem.InchIKey = inchikey;
					pubchem.InchI = inchi;
					pubchem.SMILE = smile;
					pubchem.IUPAC_trad = iupac;
					pubchem.MONOISOTOPICMASS = monomass;
					*/

					double mass = new Double(monomass);
					
					String mass_category = new Integer((int)mass).toString();
					String formula_name = elementOnly(formula);

					//appendText(str, "Indexed/" + mass_category + "/" + formula_name + ".txt");
					appendText(str, outputFolder + "/" + mass_category + "/" + formula_name + ".txt");

					if (count % 10000 == 0) {
						System.out.println(count);
					}
					count++;
				//}
			}
			in.close();
			
			System.out.println(complete.size());

		} catch (Exception e) {

		}
	}
	public static void appendText(String str, String fileName) {
		try {
		    PrintWriter out = new PrintWriter(new BufferedWriter(new FileWriter(fileName, true)));
		    out.println(str);
		    out.close();
		} catch (IOException e) {
		    //exception handling left as an exercise for the reader
		}
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

	public static LinkedList file2List(String fileName) {
		LinkedList list = new LinkedList();
		try {

			FileInputStream fstream = new FileInputStream(fileName);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				list.add(str.trim());

			}
			in.close();

		} catch (Exception e) {
			e.printStackTrace();
		}
		return list;
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

	}
}
