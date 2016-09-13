import DIGESTION.CombineFragmentLengthFile;
import DIGESTION.FragmentLengthGenerator;
import DIGESTION.PeptideDigestionSimulationLysC;
import DIGESTION.PeptideSimulatedDigestion;
import ISOTOPEDISTRIBUTION.IsotopeCalculator;
import ISOTOPEDISTRIBUTION.MS1.CalculateDegradedPeptideIsotopePattern;
import ISOTOPEDISTRIBUTION.MS1.GenerateStatisticalMatrix;
import ISOTOPEDISTRIBUTION.MS1.GenerateTheoreticalDatabase;
import METABOLOMICS.AIM.AIMFragmenter;
import METABOLOMICS.AIM.AIMFragmenterWrapper;
import METABOLOMICS.AIM.AIMPuritySampler;
import METABOLOMICS.AIM.QueryMassWithPPM;
import METABOLOMICS.AIM.QueryMassWithPPMFilter;
import METABOLOMICS.AIM.QueryStructureDatabase;
import METABOLOMIC_DATABASE.CreatePubchemDecoyDatabaseH;
import METABOLOMIC_DATABASE.CreatePubchemDecoyDatabaseNH2;
import METABOLOMIC_DATABASE.GenerateHMDBDatabase;
import METABOLOMIC_DATABASE.PubchemCombineFilesTogetherDecoyH;
import METABOLOMIC_DATABASE.PubchemCombineFilesTogetherDecoyNH2;
import METABOLOMIC_DATABASE.PUBCHEM.AppendRDBERule;
import METABOLOMIC_DATABASE.PUBCHEM.CompareTargetDecoyNumber;
import METABOLOMIC_DATABASE.PUBCHEM.DownloadPubchemXML;
import METABOLOMIC_DATABASE.PUBCHEM.GeneratePubchemIndexed;
import METABOLOMIC_DATABASE.PUBCHEM.ProcessPubChemXML;
import METABOLOMIC_DATABASE.PUBCHEM.SeparateFileListGenerateScript;
import MZXMLParser.MzXMLParser;
import MassDatabaseGeneration.AppendDatabaseInfo;
import MassDatabaseGeneration.CreateIndexFormula;
import MassDatabaseGeneration.CreateIndexFormulaExternalDB;
import MassDatabaseGeneration.CreateIndexFormulaExternalDBScript;
import MassDatabaseGeneration.CreateIndexFormulaExternalDB_H;
import MassDatabaseGeneration.CreateIndexFormulaExternalDB_NH2;
import MassDatabaseGeneration.CreateIndexFormulaOriginal;
import MassDatabaseGeneration.CreateIndexFormulaOriginalScript;
import MassDatabaseGeneration.CreateIndexFormulaScript;
import MassDatabaseGeneration.FilterFClBr;
import MassDatabaseGeneration.MassDatabaseScriptGeneration;
import PUBCHEM_STRUCTUREDATABASE.AppendPubChemTable;
import PUBCHEM_STRUCTUREDATABASE.GeneratePubchemStructureDatabase;
import PUBCHEM_STRUCTUREDATABASE.PUBCHEMStructureAddAnnotation;
import PUBCHEM_STRUCTUREDATABASE.SeparateBasedOnFormulaMass;
import PUBCHEM_STRUCTUREDATABASE.SeparateBasedOnFormulaMassFile;
import PUBCHEM_STRUCTUREDATABASE.SeparateFileBasedOnMass;
import PurityEstimation.IsotopePurityBatch;
import PurityEstimation.IsotopePurityMeasurement;
import PurityEstimation.RoughIsotopePurityMeasurement;
import PurityEstimation.RoughIsotopePurityMeasurementDTA;
import SIMULATE_FORMULA.ReconstructFormulaDBEnsureUniqFormula;
import SIMULATE_FORMULA.SimulatedFormulaPlot;
import Specialized.CreateDecoyDatabase;
import Specialized.DeIsotopeReport;
import Specialized.GeneratePubchemDecoyDatabase;
import Specialized.HMDBGenerateFragmentationListScript;
import Specialized.JUMPfIsotopeDistribution;
import Specialized.XushengDrewPaper.AppendData2FragmentCluster;
import Specialized.XushengDrewPaper.GroupCentralFormula;
import Specialized.XushengDrewPaper.MISSILEGenerateFragmentationListScript;
import Specialized.XushengDrewPaper.PaperGenerateFragmentationListScript;

/**
 * This Package is a collection of proteomics and mass spectrometry tools
 * combined as a java library
 * 
 * @author tshaw
 */
public class JumpMain {

	public static void main(String[] args) {
		try {

			if (args.length <= 0) {
				System.out.println("Not enough argument");
				printJumpJarProgramsInfo();
				System.exit(0);
			}

			String type = args[0];
			if (type.equals("-queryMassWithPPM")) {
				// System.out.println("Running Query Mass With PPM");
				/*double queryMass = new Double(args[0]);
			double tolerance = new Double(args[1]);
			String folder = args[2];
			String rgdb = args[3];
			String term = "";
			if (args.length > 4) {
				term = args[4];
			}*/
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -queryMassWithPPM [queryMass] [tolerance] [folder] [rgdb flag] [dbname (optional)]");
					System.exit(0);
				}
				runQueryMassWithPPM(args_remain);
				//
			} else if (type.equals("-QueryMassWithPPMFilter")) {
				// System.out.println("Running Isotope Calculator");
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -QueryMassWithPPMFilter [queryMass] [tolerance] [folder] [rgdb flag] [dbname (optional)]");
					System.exit(0);
				}
				QueryMassWithPPMFilter.execute(args_remain);
			} else if (type.equals("-IsotopeCalculator")) {
				// System.out.println("Running Isotope Calculator");
				String[] args_remain = getRemaining(args);
				runIsotopeCalculator(args_remain);
			} else if (type.equals("-PuritySampler")) {
				// System.out.println("Running Isotope Calculator");
				String[] args_remain = getRemaining(args);
				runAIMPuritySampler(args_remain);
			} else if (type.equals("-QueryStructureDatabase")) {
				String[] args_remain = getRemaining(args);
				QueryStructureDatabase qs = new QueryStructureDatabase();
				qs.execute(args_remain);
				//QueryStructureDatabase.execute(args_remain);
				//runQueryStructureDatabase(args_remain);
			} else if (type.equals("-FragmentSMILE")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out.println("drppm -FragmentSMILE " + AIMFragmenter.parameter_info());
					System.exit(0);
				}
				AIMFragmenter.execute(args_remain);
				// runAIMFragmenter(args_remain);
				//AIMFragmenterWrapper
			} else if (type.equals("-AIMFragmenterWrapper")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out.println("drppm -AIMFragmenterWrapper " + AIMFragmenterWrapper.parameter_info());
					System.exit(0);
				}
				AIMFragmenterWrapper.execute(args_remain);
				// runAIMFragmenter(args_remain);
				//
			} else if (type.equals("-PeptideDigestion")) {
				String[] args_remain = getRemaining(args);
				PeptideSimulatedDigestion.execute(args_remain);
			} else if (type.equals("-GenerateDatabase")) {
				String[] args_remain = getRemaining(args);
				GenerateTheoreticalDatabase.execute(args_remain);
			} else if (type.equals("-databaseStats")) {
				String[] args_remain = getRemaining(args);
				GenerateStatisticalMatrix.execute(args_remain);
			} else if (type.equals("-CreateDecoyDatabase")) {
				String[] args_remain = getRemaining(args);
				CreateDecoyDatabase.execute(args_remain);
			} else if (type.equals("-CreatePubchemDecoyDatabaseH")) {
				String[] args_remain = getRemaining(args);
				CreatePubchemDecoyDatabaseH.execute(args_remain);
			} else if (type.equals("-CreatePubchemDecoyDatabaseNH2")) {
				String[] args_remain = getRemaining(args);
				CreatePubchemDecoyDatabaseNH2.execute(args_remain);
			} else if (type.equals("-PubchemCombineFilesTogetherH")) {
				String[] args_remain = getRemaining(args);
				PubchemCombineFilesTogetherDecoyH.execute(args_remain);
			} else if (type.equals("-PubchemCombineFilesTogetherNH2")) {
				String[] args_remain = getRemaining(args);
				PubchemCombineFilesTogetherDecoyNH2.execute(args_remain);
			} else if (type.equals("-GenerateMassFormulaPubchemDecoy")) {
				String[] args_remain = getRemaining(args);
				GeneratePubchemDecoyDatabase.execute(args_remain);
			} else if (type.equals("-TrypticIsotopeGeneration")) {
				String[] args_remain = getRemaining(args);
				CalculateDegradedPeptideIsotopePattern.execute(args_remain);
			} else if (type.equals("-JUMPfIsotope")) {
				String[] args_remain = getRemaining(args);
				JUMPfIsotopeDistribution.execute(args_remain);
			} else if (type.equals("-JUMPfIsotopeReport")) {
				String[] args_remain = getRemaining(args);
				DeIsotopeReport.execute(args_remain);
			} else if (type.equals("-DownloadPubchemXML")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -DownloadPubchemXML [outputDirectory] [outputFile]");
					System.exit(0);
				}
				DownloadPubchemXML.execute(args_remain);
			} else if (type.equals("-GeneratePubchemIndexed")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -DownloadPubchemXML [outputDirectory] [outputFile]");
					System.exit(0);
				}
				GeneratePubchemIndexed.execute(args_remain);
			} else if (type.equals("-SeparateBasedOnFormulaMass")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -SeparateBasedOnFormulaMass [inputListFile] [OutputIndexFolder]");
					System.exit(0);
				}
				SeparateBasedOnFormulaMass.execute(args_remain);
				//SeparateBasedOnFormulaMassSingle
			} else if (type.equals("-SeparateBasedOnFormulaMassSingle")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -SeparateBasedOnFormulaMassSingle [inputListFile] [OutputIndexFolder]");
					System.exit(0);
				}
				SeparateBasedOnFormulaMassFile.execute(args_remain);
				//SeparateBasedOnFormulaMassFile
			} else if (type.equals("-SeparateBasedOnFormulaMassFile")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -SeparateBasedOnFormulaMassFile [inputListFile] [OutputIndexFolder]");
					System.exit(0);
				}
				SeparateBasedOnFormulaMassFile.execute(args_remain);
				//
			} else if (type.equals("-SeparateFileBasedOnMass")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -SeparateFileBasedOnMass [inputListFile] [mass divide] [OutputIndexFolder]");
					System.exit(0);
				}
				SeparateFileBasedOnMass.execute(args_remain);
				//
			} else if (type.equals("-MISSILEGenerateFragmentationListScript")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -MISSILEGenerateFragmentationListScript [inputFile]");
					System.exit(0);
				}
				MISSILEGenerateFragmentationListScript.execute(args_remain);
				// HMDBGenerateFragmentationListScript
			} else if (type.equals("-HMDBGenerateFragmentationListScript")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -HMDBGenerateFragmentationListScript [inputFile]");
					System.exit(0);
				}
				HMDBGenerateFragmentationListScript.execute(args_remain);
				// PaperGenerateFragmentationListScript
			} else if (type.equals("-PaperGenerateFragmentationListScript")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -PaperGenerateFragmentationListScript [inputFile]");
					System.exit(0);
				}
				PaperGenerateFragmentationListScript.execute(args_remain);
				// AppendData2FragmentCluster
			} else if (type.equals("-AppendData2FragmentCluster")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -AppendData2FragmentCluster [inputFile] [fragmentFile] [hmdbFile] [outputFile]");
					System.exit(0);
				}
				AppendData2FragmentCluster.execute(args_remain);
				//
			} else if (type.equals("-CompareTargetDecoy")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -CompareTargetDecoy [inputFile] [outputFile]");
					System.exit(0);
				}
				CompareTargetDecoyNumber.execute(args_remain);
				//
			} else if (type.equals("-ReconstructFormulaDBEnsureUniqFormula")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -ReconstructFormulaDBEnsureUniqFormula [inputFolder] [outputFolder]");
					System.exit(0);
				}
				ReconstructFormulaDBEnsureUniqFormula.execute(args_remain);
				// SimulatedFormulaPlot
			} else if (type.equals("-SimulatedFormulaPlot")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -SimulatedFormulaPlot [inputFolder] [outputFile]");
					System.exit(0);
				}
				SimulatedFormulaPlot.execute(args_remain);
				// SimulatedFormulaPlot
			} else if (type.equals("-MassDatabaseScriptGeneration")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -MassDatabaseScriptGeneration [output_mass_prefix] [output_shell_prefix] [num of divided tasks] [formulaType]");
					System.exit(0);
				}
				MassDatabaseScriptGeneration.execute(args_remain);
				// CreateIndexFormulaScript
			} else if (type.equals("-CreateIndexFormulaScript")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -CreateIndexFormulaScript [inputFile] [outputFolder] [min] [max] [freq] [outputFile]");
					System.exit(0);
				}
				CreateIndexFormulaScript.execute(args_remain);
				//
			} else if (type.equals("-CreateIndexFormula")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -CreateIndexFormula [inputFile] [outputFolder] [min] [max]");
					System.exit(0);
				}
				CreateIndexFormula.execute(args_remain);

			} else if (type.equals("-CreateIndexFormulaOriginal")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -CreateIndexFormulaOriginal [inputFile] [outputFolder] [min] [max]");
					System.exit(0);
				}
				CreateIndexFormulaOriginal.execute(args_remain);
				// CreateIndexFormulaOriginalScript
			} else if (type.equals("-CreateIndexFormulaOriginalScript")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -CreateIndexFormulaOriginalScript [inputFile] [outputFolder] [parallel_num]");
					System.exit(0);
				}
				CreateIndexFormulaOriginalScript.execute(args_remain);
				// CreateIndexFormulaExternalDB
			} else if (type.equals("-CreateIndexFormulaExternalDBScript")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {

					System.out
							.println("jumpjar -CreateIndexFormulaExternalDBScript [inputFile] [outputFolder] [parallel_num] [externalFiles] [databaseNames] [databaseIndex] [NormalDecoyFlag] [ValenceRule]");
					System.exit(0);
				}
				CreateIndexFormulaExternalDBScript.execute(args_remain);
				//
			} else if (type.equals("-CreateIndexFormulaExternalDB")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {

					System.out
							.println("jumpjar -CreateIndexFormulaExternalDB [inputFile] [folder] [min] [max] [otherReference] [referenceNames] [formulaColIndex] [NormalDecoyFlag] [ValenceRule]");
					System.exit(0);
				}
				CreateIndexFormulaExternalDB.execute(args_remain);
				// CreateIndexFormulaExternalDB_H
			} else if (type.equals("-CreateIndexFormulaExternalDB_H")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {

					System.out
							.println("jumpjar -CreateIndexFormulaExternalDB_H [inputFile] [folder] [min] [max] [otherReference] [referenceNames] [formulaColIndex] [NormalDecoyFlag] [ValenceRule]");
					System.exit(0);
				}
				CreateIndexFormulaExternalDB_H.execute(args_remain);
				// CreateIndexFormulaExternalDB_H
			} else if (type.equals("-CreateIndexFormulaExternalDB_NH2")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {

					System.out
							.println("jumpjar -CreateIndexFormulaExternalDB_NH2 [inputFile] [folder] [min] [max] [otherReference] [referenceNames] [formulaColIndex] [NormalDecoyFlag] [ValenceRule]");
					System.exit(0);
				}
				CreateIndexFormulaExternalDB_NH2.execute(args_remain);
				// CreateIndexFormulaExternalDB_H
			} else if (type.equals("-FilterFClBr")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {

					System.out
							.println("jumpjar -FilterFClBr " + FilterFClBr.parameter_info());
					System.exit(0);
				}
				FilterFClBr.execute(args_remain);
				//
			} else if (type.equals("-AppendDatabaseInfo")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					/*
					 * String referenceFile = args[0]; int massColIndex = new
					 * Integer(args[1]); String databasePath = args[2]; String
					 * databaseName = args[3]; String test_flag_str = args[4];
					 */
					System.out
							.println("jumpjar -AppendDatabaseInfo [referenceFile] [referenceFormulaIndex] [databasePath] [databaseName] [testFlag]");
					System.exit(0);
				}
				AppendDatabaseInfo.execute(args_remain);
				//
			} else if (type.equals("-SeparateFileList")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -SeparateFileList [inputFile] [num]");
					System.exit(0);
				}
				SeparateFileListGenerateScript.execute(args_remain);
				//
			} else if (type.equals("-GroupCentralFormula")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -GroupCentralFormula [inputFile] [matrixFile]");
					System.exit(0);
				}
				GroupCentralFormula.execute(args_remain);
				//GenerateHMDBDatabase
			} else if (type.equals("-GenerateHMDBDatabase")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -GenerateHMDBDatabase " + GenerateHMDBDatabase.parameter_info());
					System.exit(0);
				}
				GenerateHMDBDatabase.execute(args_remain);
				//PUBCHEMStructureAddAnnotation
			} else if (type.equals("-PUBCHEMStructureAddAnnotation")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -PUBCHEMStructureAddAnnotation " + PUBCHEMStructureAddAnnotation.parameter_info());
					System.exit(0);
				}
				PUBCHEMStructureAddAnnotation.execute(args_remain);
				//PeptideDigestionSimulationLysC
			} else if (type.equals("-PeptideDigestionSimulationLysC")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -PeptideDigestionSimulationLysC " + PeptideDigestionSimulationLysC.parameter_info());
					System.exit(0);
				}
				PeptideDigestionSimulationLysC.execute(args_remain);
				//FragmentLengthGenerator
			} else if (type.equals("-FragmentLengthGenerator")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -FragmentLengthGenerator " + FragmentLengthGenerator.parameter_info());
					System.exit(0);
				}
				FragmentLengthGenerator.execute(args_remain);
				//CombineFragmentLengthFile
			} else if (type.equals("-CombineFragmentLengthFile")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -CombineFragmentLengthFile " + CombineFragmentLengthFile.parameter_info());
					System.exit(0);
				}
				CombineFragmentLengthFile.execute(args_remain);
				// ProcessPubChemXML
			} else if (type.equals("-ProcessPubChemXML")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -ProcessPubChemXML " + ProcessPubChemXML.parameter_info());
					System.exit(0);
				}
				ProcessPubChemXML.execute(args_remain);
				// AppendRDBERule
			} else if (type.equals("-AppendRDBERule")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -AppendRDBERule " + AppendRDBERule.parameter_info());
					System.exit(0);
				}
				AppendRDBERule.execute(args_remain);
				// AppendPubChemTable
			} else if (type.equals("-AppendPubChemTable")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -AppendPubChemTable " + AppendPubChemTable.parameter_info());
					System.exit(0);
				}
				AppendPubChemTable.execute(args_remain);
				// GeneratePubchemStructureDatabase
			} else if (type.equals("-GeneratePubchemStructureDatabase")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -GeneratePubchemStructureDatabase " + GeneratePubchemStructureDatabase.parameter_info());
					System.exit(0);
				}
				GeneratePubchemStructureDatabase.execute(args_remain);
				// MzXMLParser
			} else if (type.equals("-MzXMLParser")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -MzXMLParser " + MzXMLParser.parameter_info());
					System.exit(0);
				}
				MzXMLParser.execute(args_remain);
				// IsotopePurityMeasurement
			} else if (type.equals("-IsotopePurityMeasurement")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -IsotopePurityMeasurement " + IsotopePurityMeasurement.parameter_info());
					System.exit(0);
				}
				IsotopePurityMeasurement.execute(args_remain);
				// IsotopePurityBatch
			} else if (type.equals("-IsotopePurityBatch")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -IsotopePurityBatch " + IsotopePurityBatch.parameter_info());
					System.exit(0);
				}
				IsotopePurityBatch.execute(args_remain);
				// RoughIsotopePurityMeasurement
			} else if (type.equals("-RoughIsotopePurityMeasurement")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -RoughIsotopePurityMeasurement " + RoughIsotopePurityMeasurement.parameter_info());
					System.exit(0);
				}
				RoughIsotopePurityMeasurement.execute(args_remain);
				// RoughIsotopePurityMeasurementDTA
			} else if (type.equals("-RoughIsotopePurityMeasurementDTA")) {
				String[] args_remain = getRemaining(args);
				if (args_remain.length == 0) {
					System.out
							.println("jumpjar -RoughIsotopePurityMeasurementDTA " + RoughIsotopePurityMeasurementDTA.parameter_info());
					System.exit(0);
				}
				RoughIsotopePurityMeasurementDTA.execute(args_remain);
				// 
			} else {
				System.out.println("Here are the programs available");
				printJumpJarProgramsInfo();
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static void printJumpJarProgramsInfo() {

		System.out
				.println("-queryMassWithPPM a program for identifying formulas based on mass");
		System.out
				.println("-IsotopeCalculator a program for calculating the isotope given a formula");
		System.out.println("-PuritySampler");
		System.out.println("-QueryStructureDatabase");
		System.out.println("-FragmentSMILE");
		System.out.println("-PeptideDigestion");
		System.out.println("-GenerateDatabase");
		System.out
				.println("-CompareTargetDecoy a program for comparing target and decoy equations");
		System.out
				.println("-ReconstructFormulaDBEnsureUniqFormula a program for ensuring uniq formula in FormulaDB");
		System.out
				.println("-SimulatedFormulaPlot generate a matrix file with information about the database");
	}

	public static String[] getRemaining(String[] args) {
		String[] args_remain = new String[args.length - 1];
		for (int i = 1; i < args.length; i++) {
			args_remain[i - 1] = args[i];
		}
		return args_remain;
	}

	public static void runQueryMassWithPPM(String[] args) {
		QueryMassWithPPM qm = new QueryMassWithPPM();
		qm.execute(args);
	}

	public static void runQueryStructureDatabase(String[] args) {
		QueryStructureDatabase qs = new QueryStructureDatabase();
		qs.execute(args);
	}

	public static void runIsotopeCalculator(String[] args) {
		IsotopeCalculator ic = new IsotopeCalculator();
		ic.execute(args);
	}

	public static void runAIMPuritySampler(String[] args) {
		AIMPuritySampler aim_sampler = new AIMPuritySampler();
		aim_sampler.execute(args);
	}
	/*
	 * public static void runAIMFragmenter(String[] args) { AIMFragmenter
	 * aim_fragmenter = new AIMFragmenter();
	 * 
	 * aim_fragmenter.execute(args); }
	 */
}
