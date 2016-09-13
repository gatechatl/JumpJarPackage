package Specialized;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class HMDBGenerateFragmentationListScript {

	public static void execute(String[] args) {
		try {
			String inputFile = args[0];
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				if (split.length >= 5) {
					String id = split[0];
					String smile = split[4];
					smile = smile.replaceAll("\\(", "\\\\(");
					smile = smile.replaceAll("\\)", "\\\\)");
					System.out.println("perl JumpJar.pl -program=SmileFragmenter -Smile='" + smile + "' -mass=50 -depth=2 > HMDB/fragment_result_" + id);
				} else {
					System.out.println("BAD: " + str);
				}
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
