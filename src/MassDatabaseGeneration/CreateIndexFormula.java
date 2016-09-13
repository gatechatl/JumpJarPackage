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


/**
 * Convert the files into load database script
 * @author tshaw
 *
 */
public class CreateIndexFormula {

	public static void execute(String[] args) {
		
		try {
			HashMap map = new HashMap();
			String folderName = args[0];
			String outputFolder = args[1];
			File folder = new File(folderName);
			for (File file: folder.listFiles()) {
				String inputFile = file.getPath();
				String name = inputFile.replaceAll(".txt", "");
				map.put(name, name);
				
				String outputFile = outputFolder + "/" + inputFile + ".dat";
				FileWriter fwriter = new FileWriter(outputFile);
				BufferedWriter out = new BufferedWriter(fwriter);
				
				FileInputStream fstream = new FileInputStream(inputFile);
				DataInputStream din = new DataInputStream(fstream);
				BufferedReader in = new BufferedReader(new InputStreamReader(din));
				while (in.ready()) {
					String str = in.readLine();
					String[] split = str.split(":");
					out.write("\"" + split[0] + "\",\"" + split[1] + "\"\n");
					// get this in the right format for each file
				}
				in.close();				
				out.close();
			}
			
			String outputFile = "index_names.txt";
			FileWriter fwriter = new FileWriter(outputFile);
			BufferedWriter out = new BufferedWriter(fwriter);
			
			Iterator itr = map.keySet().iterator();
			while (itr.hasNext()) {
				String key = (String)itr.next();
				out.write(key + "\n");
			}
			out.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}

