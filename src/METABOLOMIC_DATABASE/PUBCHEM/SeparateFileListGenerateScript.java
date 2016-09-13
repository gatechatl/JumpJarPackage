package METABOLOMIC_DATABASE.PUBCHEM;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;

/**
 * To Assist with parallelization
 * We separate the pubchem FileList into 
 * @author tshaw
 *
 */
public class SeparateFileListGenerateScript {

	public static void execute(String[] args) {
		
		try {
			String inputFile = args[0];
			int num = new Integer(args[1]);
			String outputFolder = args[2];
			String outputScriptFile = args[3];
			
			FileWriter fwriter = new FileWriter(outputScriptFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			File file = new File(outputFolder);
			if (!file.exists()) {
				file.mkdir();
			}
			for (int i = 1; i <= num; i++) {
				file = new File(inputFile + "_" + i);
				out.write("jumpjar -SeparateBasedOnFormulaMass " + inputFile + "_" + i + " " + outputFolder + "\n");
				if (file.exists()) {
					file.delete();
				}
			}
			out.close();
						
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {				
				for (int i = 1; i <= num; i++) {
					if (in.ready()) {
						String str = in.readLine();
						appendText(str, inputFile + "_" + i);
					}
				}			
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
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
}
