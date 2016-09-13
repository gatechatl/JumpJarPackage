package MassDatabaseGeneration;

public class CreateIndexFormulaExternalDBScript {
	public static void execute(String[] args) {
		
		String inputFile = args[0];
		String outputFolder = args[1];
		int parallel_num = new Integer(args[2]);
		String externalFiles = args[3];
		String databaseNames = args[4];
		String databaseIndex = args[5];
		String generateDecoy = args[6]; // N or D flag
		String valenceRule = args[7];
		int buffer = 1600 / parallel_num;
		for (int i = 0; i <= 1600; i = i + buffer) {
			int max = i + buffer;
			if (max > 1600) {
				max = 1600;
			}
			System.out.println("jumpjar -CreateIndexFormulaExternalDB_H " + inputFile + " " + outputFolder + " " + i + " " + (i + buffer) + " " + externalFiles + " " + databaseNames + " " + databaseIndex + " " + generateDecoy + " " + valenceRule);
		}
	}
}
