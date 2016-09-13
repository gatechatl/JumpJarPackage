package Specialized.XushengDrewPaper;

import java.io.BufferedReader;
import java.io.DataInputStream;
import java.io.FileInputStream;
import java.io.InputStreamReader;

public class PaperGenerateFragmentationListScript {
	public static void execute(String[] args) {
		try {
			String inputFile = args[0];
			int index = 1;
			FileInputStream fstream = new FileInputStream(inputFile);
			DataInputStream din = new DataInputStream(fstream);
			BufferedReader in = new BufferedReader(new InputStreamReader(din));
			while (in.ready()) {
				String str = in.readLine();
				String[] split = str.split("\t");
				String id = "ID" + index;
				String smile = split[1];
				//smile = smile.replaceAll("\\(", "\\\\(");
				//smile = smile.replaceAll("\\)", "\\\\)");
				System.out.println("perl JumpJar.pl -program=SmileFragmenter -Smile='" + smile + "' -mass=50 -depth=2 > fragment_result_" + id);
				index++;
			}
			in.close();
			
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
