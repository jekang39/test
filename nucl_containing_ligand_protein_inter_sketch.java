
import java.io.*;
import java.util.*;
import java.util.regex.*;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.FileUtils;
import org.apache.commons.lang3.ArrayUtils;
import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;

public class  nucl_containing_ligand_protein_inter_sketch {
    
public static void main(String[] args) throws IOException,InterruptedException{

    nucl_containing_ligand_protein_inter_sketch aa = new nucl_containing_ligand_protein_inter_sketch();
    String prot_domain_binding_profile = "ADEP";
    String dna_binding_site = "DYCACCGAC";
    aa.get_binding_profiles_human(prot_domain_binding_profile,dna_binding_site);
    //instantiate inner class NBDB interface  
    nucl_containing_ligand_protein_inter_sketch.NBDB nbdb = new nucl_containing_ligand_protein_inter_sketch.NBDB();
    ////inner class functions
    nbdb.make_uniq_ligand_list();
    nbdb.make_profile_list();
    nbdb.make_ligand_name_list();
    nbdb.map_pos_freq_matrix();
    nbdb.map_profile_ligand_name();
}


//LinkedHashMap<String,LinkedList<String>> tf_to_dna_binding_domain_human;
LinkedHashMap<String,LinkedHashMap<String,LinkedHashMap<String,String>>> tf_to_domain_binding_motif_human;
LinkedHashMap<String,LinkedHashMap<String,LinkedList<String>>> tf_to_domain_binding_motif_matrix_human;

//key1:cis-bp,key2:domain-seq,key3:motif-id
int aa_group_trimer_array_length;
int dna_group_trimer_array_length;
public void  get_binding_profiles_human(String domain_seq, String motif_seq)  throws IOException{

String[][] aa_group_trimer_array= {{"PPP","PPP","PPP","PPP","PPP","PPP"},
{"SSS","SSS","SSS","SSS","SSS","SSS"},
{"NNN","NNN","NNN","NNN","NNN","NNN"},
{"HHH","HHH","HHH","HHH","HHH","HHH"},
{"RRR","RRR","RRR","RRR","RRR","RRR"},
{"PPN","PNP","PPN","PNP","NPP","NPP"},
{"PNN","PNN","NPN","NNP","NPN","NNP"},
{"PPS","PSP","PPS","PSP","SPP","SPP"},
{"PSS","PSS","SPS","SSP","SPS","SSP"},
{"PPH","PHP","PPH","PHP","HPP","HPP"},
{"PHH","PHH","HPH","HHP","HPH","HHP"},
{"PPR","PRP","PPR","PRP","RPP","RPP"},
{"PRR","PRR","RPR","RRP","RPR","RRP"},
{"PNS","PSN","NPS","NSP","SPN","SNP"},
{"PNH","PHN","NPH","NHP","HPN","HNP"},
{"PNR","PRN","NPR","NRP","RPN","RNP"},
{"PSH","PHS","SPH","SHP","HPS","HSP"},
{"PSR","PRS","SPR","SRP","RPS","RSP"},
{"PHR","PRH","HPR","HRP","RPH","RHP"},
{"NNS","NSN","NNS","NSN","SNN","SNN"},
{"NNH","NHN","NNH","NHN","HNN","HNN"},
{"NSS","NSS","SNS","SSN","SNS","SSN"},
{"NSH","NHS","SNH","SHN","HNS","HSN"},
{"NSR","NRS","SNR","SRN","RNS","RSN"},
{"NHR","NRH","HNR","HRN","RNH","RHN"},
{"NRN","NNR","RNN","RNN","NNR","NRN"},
{"NRR","NRR","RNR","RRN","RNR","RRN"},
{"SSH","SHS","SSH","SHS","HSS","HSS"},
{"SSR","SRS","SSR","SRS","RSS","RSS"},
{"SHH","SHH","HSH","HHS","HSH","HHS"},
{"SHR","SRH","HSR","HRS","RSH","RHS"},
{"SRR","SRR","RSR","RRS","RSR","RRS"},
{"HNH","HHN","NHH","NHH","HHN","HNH"},
{"HHR","HRH","HHR","HRH","RHH","RHH"},
{"HRR","HRR","RHR","RRH","RHR","RRH"}
};

String[][] dna_group = {
 {"A","C","G","T"},/*0*/
 {"M","R","W","S","Y","K"},/*1*/
 {"V","H","D","B"},/*2*/
 {"N"}/*3*/,
 
 };
 
 int[][] aa_perm = {{0,1,2}, {0,2,1}, {1,0,2}, {1,2,0}, {2,0,1}, {2,1,0}};
 
 int[][] permutations={ {0,1,2}, {0,2,1}, {0,1,3}, {0,3,1}, {0,2,3}, {0,3,2}, {1,0,2}, {1,2,0}, {1,0,3}, {1,3,0}, {1,2,3}, {1,3,2}, {2,0,1}, {2,1,0}, {2,0,3}, {2,3,0}, {2,1,3}, {2,3,1}, {3,0,1}, {3,1,0}, {3,0,2}, {3,2,0}, {3,1,2}, {3,2,1}};
 
 String[] dna_group_letter = {"Z","O","T","N"};
 
 String[] dna_one_letter = {"A"	,"C"	,"G","T","M","R","W","S","Y","K"	,"V"	,"H"	,"D"	,"B"	,"N"};
 
 int[] dna_one_letter_group_int = {0,0,0,0,1,1,1,1,1,1,2,2,2,2,3};
 
 	String[][] aa_map ={{"Ala","A"},{"Arg","R"},{"Asn","N"},{"Asp","D"},{"Cys","C"},{"Gln","Q"},{"Glu","E"},{"Gly","G"},{"His","H"},{"Ile","I"},{"Leu","L"},{"Lys","K"},{"Met","M"},{"Phe","F"},{"Pro","P"},{"Ser","S"},{"Thr","T"},{"Trp","W"},{"Tyr","Y"},{"Val","V"}}; 
 	
 	String[][] aa_group ={{"R","K","S","T"}/*P:POSITIVE, POLAR-UNCHARGED 0*/,
{"D","E","N","Q"}/*N:NEGATIVE;POLAR-UNCHARGED 1*/,
{"C","G","P","H"}/*S:SPECIAL 2*/,
{"A","V","I","L","M"}/*H:HYDROPHBIC 3*/,
{"F","W","Y"}/*R:RING 4*/
}; 

 	/*it will increase sensitivity if we convert each special aa into a separate group
String[][] aa_group ={
{"R","K","S","T"},
{"D","E","N","Q"},

{"A","V","I","L","M"},
{"F","W","Y"},
{"C"},
{"G"}
{"P"}
{"H"}
}; 
*/
  	String[] aa_group_letter = {"P","N","S","H","R"};//
  	String[] three_letter = {"VAL","LEU","LEU","GLU","GLY","PRO","PRO","HIS","SER","GLY","LYS","THR","ALA","LEU","ALA","ALA","LYS","ILE","ALA","GLU","GLU","SER","ASN","PHE","PRO","PHE","ILE","LYS","ILE","CYS","SER","PRO","ASP","LYS","MET","ILE","GLY","PHE","SER","GLU","THR","ALA","LYS","CYS","GLN","ALA","MET","LYS","LYS","ILE","PHE","ASP","ASP","ALA","TYR","LYS","SER","GLN","LEU","SER","CYS","VAL","VAL","VAL","ASP","ASP","ILE","GLU","ARG","LEU","LEU","ASP","TYR","VAL","PRO","ILE","GLY","PRO","ARG","PHE","SER","ASN","LEU","VAL","LEU","GLN","ALA","LEU","LEU","VAL","LEU","LEU","LYS","LYS","ALA","PRO","PRO","GLN","GLY","ARG","LYS","LEU","LEU","ILE","ILE","GLY","THR","THR","SER","ARG","LYS","ASP","VAL","LEU","GLN","GLU","MET","GLU","MET","LEU","ASN","ALA","PHE","SER","THR","THR","ILE","HIS","VAL","PRO","ASN"};
	String[] one_letter = {"A","R","N","D","C","Q","E","G","H","I","L","K","M","F","P","S","T","W","Y","V"};
	int[] one_letter_group_int = {3,0,1,1,2,1,1,2,2,3,3,0,3,4,2,0,0,4,4,3 };
	
	String[] dna_group_trimer_array = new String[64];
	dna_group_trimer_array_length = dna_group_trimer_array.length;
 	int count = 0;
 	
 	for(int a = 0; a < dna_group_letter.length; a++){
 	
 		for(int b = 0; b < dna_group_letter.length; b++){
 			
 			for(int c = 0; c < dna_group_letter.length; c++){
 			
 				dna_group_trimer_array[count] = dna_group_letter[a] + dna_group_letter[b] + dna_group_letter[c];
 				count = count +1;
 			}
 		}
 	}
 			
 	System.out.println("dna_group_trimer_array:" + count );
 	System.out.println("dna trimer:"+ Arrays.toString(dna_group_trimer_array));

	System.out.println("get_binding_profiles_human():" );
	

	String domain_matrix = "";
	String motif_matrix = "";
		
	
	aa_group_trimer_array_length = aa_group_trimer_array.length;
			
	int[] aa_abs_pres = new int[aa_group_trimer_array.length];
			
	for(int i = 0; i < domain_seq.length(); i++){
			
		if((i+3)> domain_seq.length()){
				
			break;
		}
		String trimer_temp = domain_seq.substring(i,i+3);
				
		int aa_group_ind = -1;
				
		for(int h = 0; h < aa_group_trimer_array.length; h++){
				
			for(int e = 0; e < aa_group_trimer_array[h].length; e++){
				
				String convert_group = "";
						
				for(int r = 0; r < trimer_temp.length(); r++){
						
					String temp = trimer_temp.substring(r,r+1);
						
					boolean flag = false;
							
					for(int t = 0; t < aa_group.length; t++){
							
						for(int y = 0; y < aa_group[t].length; y++){
							
							if(temp.equals(aa_group[t][y])){
									
								convert_group = convert_group + aa_group_letter[t];
								flag = true;
								break;
							}
						}
								
						if(flag){
								
							break;
						}
					}
				}
							
				System.out.println("convert_group:" + convert_group );
				if(convert_group.equals(aa_group_trimer_array[h][e])){
						
					aa_group_ind = h;
					aa_abs_pres[aa_group_ind] =aa_abs_pres[aa_group_ind]+1;
					break;
						
				}
			}//e
					
				
		}//h
		System.out.println("aa_group_ind:"+ aa_group_ind + ":" + trimer_temp);
				
	}//i
			
	domain_matrix = "";
			
	for(int t = 0; t < aa_abs_pres.length; t++){
			
		domain_matrix = domain_matrix + aa_abs_pres[t] + ",";
	}
			
	domain_matrix = domain_matrix.substring(0,domain_matrix.length()-1);
		
	//aa matrix done
			
	//dna motif matrix start
	int[] dna_abs_pres = new int[dna_group_trimer_array.length];
	for(int i = 0; i < motif_seq.length(); i++){
												
		if((i+3)> motif_seq.length()){
													
			break;
		}
		String trimer_temp = motif_seq.substring(i,i+3);
													
		int dna_group_ind = -1;
					
		for(int h = 0; h < dna_group_trimer_array.length; h++){
													
			if(trimer_temp.equals(dna_group_trimer_array[h])){
														
				dna_group_ind = h;
				dna_abs_pres[dna_group_ind] = dna_abs_pres[dna_group_ind]+1;
				break;
														
														
			}
													
		}
					
	}//i
				
				
	motif_matrix = "";
				
	for(int t = 0; t < dna_abs_pres.length; t++){
				
		motif_matrix = motif_matrix+ dna_abs_pres[t] + ",";
	}
				
	motif_matrix = motif_matrix.substring(0,motif_matrix.length()-1);
				
	//apply NBDB information to domain and motif	
	
}//method

//NBDB
static class NBDB{

	String[] ligand_type = {"RBP","RBPF","RBPN","RBPS","RBPSO","RBSO","TOP","OP","TP","RPF","RPFO"};

	String[][] ligand_name = {{"CTP","CMP", "GMP","cAMP","AMP","ATP","GDP","PAP","cGMP","ADP","GTP","c-di-GMP","c-di-AMP"},{"FAD"},{"NAD","NADP","NAP"},{"PPS"},{"CoA","Acetyl-CoA"},{"SAM"},{"THD"},{"PLP"},{"ThPP", "TPP"},{"FMN"},{"F420", "F42"}};

	String[] ligand_name_indi = {"CTP","CMP", "GMP","cAMP","AMP","ATP","GDP","PAP","cGMP","ADP","GTP","c-di-GMP","c-di-AMP","FAD","NAD","NADP","NAP","PPS","CoA","Acetyl-CoA","SAM","THD","PLP","ThPP","TPP", "FMN","F420", "F42"};

	//{{CTP,GMP,cAMP,AMP,c-di-GMP,ATP,GDP,PAP,c-di-AMP,cGMP,ADP,GTP},{FAD(H)},{NADP,NAD(H)},{PPS},{Acetyl-CoA,CoA},{SAM},{THD},{PLP},{ThPP},{FMN},{F420}};

	String[] profile = {"ADEP","ADLG","AGDY","ANNAGH","ANRGG","APGSYRD","ATG","AVGR","CDDGRT","CFY","CNCYC","CPYN","DAAN","DAHG","DAYQGF","DDQKPPG","DED","DEY","DFDHYD","DFRRNGEARLF","DGFPQ","DGGGTD","DGGIFF","DGGSQ","DGKE","DGTD","DGTS","DI","DKP","DLFDHTAF","DLVVGGR","DNGKG","DPDGS","DPYPKD","DQ","DQFP","DRAVP","DRDR","DRMEYE","DRP","DRRDFTN","DSYF","DYFC","DYPPFD","EFGP","EFGTP","EFSHR","EFWIEENWTD","EGGSG","EGHPDD","EGVEC","EKLS","EKNRA","ELKMDGIA","END","ENPDW","ENR","EPGR","ESTGFH","ETAL","ETKPDR","EYGFGEF","EYGG","FDL","FDQKRHYKRQ","FDTP","FESA","FGEAAKA","FGYDR","FMKNN","FPGGDK","FTGTE","FVWW","FYPTY","G4FG","GAQTG","GAR","GAVD","GCGSD","GDDFGTGYSS","GDFGWTNI","GDGMET","GDGTT","GDNNAG","GDRAG","GEETTTG","GESGGKT","GFEGT","GFGGGG","GGEP","GGGTGG","GGGTP","GGHNS","GGI","GGIG","GGR","GGRK","GGSGP","GGTD","GGTRP","GHGS","GKSTN","GKW","GMGPG","GNKP","GNKTD","GNPRFR","GPG","GPKPRP","GQGNQ","GRFG","GRGMR","GRGRP","GRGSNF","GRLNHM","GRP","GRTM","GTDGR","GTGSSEMKW","GTHP","GVCDVD","GVVV","GxGGxG","GxGxxG","GxxGxG","HDEE","HDN","HDRP","HGAKKVP","HGQHGR","HHHAH","HIGH","HPSFN","HPTE","HSGDR","HTFGNP","HTGRPNQNIP","HVEKP","HVVRRGYE","HWH","IANPTEGRH","KADFL","KECY","KENSKR","KFYGG","KGRG","KKE","KVVTGHEPW","KYGHG","LAGSVLG","LDGCGG","LFFGCR","LFGVGEP","LNLDNGT","LTYYM","MDHFF","MNAGCT","MTHTQD","NDDD","NDHSLF","NIGH","NKYEGPRYYG","NNGG","NPGFGGPD","NPGT","NQGDEP","NRLSD","NTKRG","NVTAFW","PDCGGGSNGF","PDKG","PFGGEDF","PGGDR","PGGTW","PGKPP","PGLH","PGPHNGKIT","PKPG","PLVGNKDL","PMNKFD","PNFGT","PNPG","PNQECH","PRHP","PSNQW","PVGI","PWFQ","PYFGGFN","QKGEPYM","QRLTRFDVR","RGTPGTIR","RKGGDW","RLFCCC","RNPDPY","RPDKQVTY","RPTRR","RQSGTGAF","RRSR","RWHTCA","SD","SDGHR","SGGDS","SGKTRG","SGLPR","SGNAG","SGRGDKD","SIAGEGLFK","SKYRG","SNHGGR","SQASQS","SRRR","TAK","TEPGSD","TFFVHKKKPM","TGEN","TGGAD","TGGSSGIG","TLGTR","TNGLL","TPFFE","TTGYGG","TTRREGYF","TYAGD","VCCKNNCHST","VDGT","VGGL","VHDW","VND","VPDY","VRRPE","VSRGKW","WDGKE","WDPDP","WDTAGQE","WFQFW","WFQY","WGGP","WHHH","WYDNE","WYPHPP","WYRP","YASK","YFFLM","YGAF","YGHGGGWG","YGLR","YHNH","YQPELRW"};

	String[] matrix_head = {"A", "C", "D", "E", "F", "G", "H", "I", "K", "L", "M", "N", "P", "Q", "R", "S", "T", "V", "W", "Y"};

	LinkedList<String> uniq_ligand_list;
	String[] uniq_ligand;
	public void make_uniq_ligand_list() throws IOException{

		uniq_ligand_list = new LinkedList();
		for(int e = 0; e < ligand_type.length; e++){

			for(int w = 0; w < ligand_type[e].length(); w++){

				String indi = ligand_type[e].substring(w,w+1);
				if(!uniq_ligand_list.contains(indi)){
					uniq_ligand_list.add(indi);
				}
			} 
			
		}
		uniq_ligand = uniq_ligand_list.toArray(new String[0]);
		System.out.println("uniq_ligand:" + Arrays.toString(uniq_ligand));
		//R, B, P, F, N, S, O, T
	}

	LinkedList<String> profile_list;
	public void make_profile_list() throws IOException{

		profile_list = new LinkedList();
		for(int e = 0; e < profile.length; e++){
			profile_list.add(profile[e]);
		}
	}

	LinkedList<String> ligand_name_list;
	public void make_ligand_name_list() throws IOException{

		ligand_name_list = new LinkedList();
		for(int e = 0; e < ligand_name_indi.length; e++){

			ligand_name_list.add(ligand_name_indi[e]);
			
		}
	}

	LinkedHashMap<String,Double[][]> profile_pos_freq_matrix;//Double[][]: position.matrix_head 
	LinkedHashMap<String,Double[][]> profile_pos_norm_matrix;

	public void map_pos_freq_matrix() throws IOException{

		profile_pos_freq_matrix = new LinkedHashMap();
		profile_pos_norm_matrix = new LinkedHashMap();

		String input="art_foundation_log_online/data/pdb_hits.txt";
		//ADEP	1DXY:A	1.61	244	LSNLKSGKLAGVGIDTYEYETEDLLNLAKH	101-299	51735	51830	c.2.1	c.2.1.4	NAD(P)-binding Rossmann-fold domains	Formate/glycerate dehydrogenases, NAD-domain	(NAD_N7N|ASP|258)
		//ADEP	1HKU:A	4.2	265	AQALKEGRIRGAALDVHESEPFSFSQGPLK	115-307	51735	51830	c.2.1	c.2.1.4	NAD(P)-binding Rossmann-fold domains	Formate/glycerate dehydrogenases, NAD-domain	(NAD_N7N|ASP|279)

		Pattern space_pattern = Pattern.compile("\\s+");
		Pattern tab_pattern = Pattern.compile("\\t+");
		Pattern comma_pattern = Pattern.compile(",");
		Pattern colon_pattern = Pattern.compile(":");

		List<String> lines  = FileUtils.readLines(new File(input));
		
		for(int q= 0; q < lines.size(); q++){
		
			String one_line = lines.get(q).trim();
			
			String[] space_split = space_pattern.split(one_line);
			int prof_ind = -1;
			String prof_name = "";
			if(profile_list.contains(space_split[0])){

				for(int e = 0; e < profile.length; e++){
					if(profile[e].equals(space_split[0])){
						prof_ind = e;
						prof_name = space_split[0];
						break;
					}
				}
			}

			if(prof_ind != -1){

				if(profile_pos_freq_matrix.containsKey(prof_name)){

					Double[][] exist = profile_pos_freq_matrix.get(prof_name);

					for(int x = 0; x < space_split[4].length(); x++){

						String indi = space_split[4].substring(x,x+1);

						for(int y = 0; y < matrix_head.length; y++){

							if(indi.equals(matrix_head[y])){

								exist[x][y] = exist[x][y]  +new Double(1.0);
								break;
							}
						}
					}
					profile_pos_freq_matrix.put(prof_name,exist);

				}else{
		
					Double[][] exist = new Double[space_split[4].length()][matrix_head.length];
					for(int x = 0; x < space_split[4].length(); x++){

						for(int y = 0; y < matrix_head.length; y++){
							exist[x][y] = new Double(0.0);
						}
					}

					for(int x = 0; x < space_split[4].length(); x++){

						String indi = space_split[4].substring(x,x+1);

						for(int y = 0; y < matrix_head.length; y++){

							if(indi.equals(matrix_head[y])){

								exist[x][y] = new Double(1.0);
								break;
							}
						}
					}
					profile_pos_freq_matrix.put(prof_name,exist);
				}//else
			}//!=-1
		}//q

		List<String> keys_list = new LinkedList<String>(profile_pos_freq_matrix.keySet());	

		for(int q = 0; q < keys_list.size();	q++){

			Double[][] exist = profile_pos_freq_matrix.get(keys_list.get(q));
			Double[][] norm_obj = new Double[exist.length][matrix_head.length];
			double all_max = 0;
			double all_min = 0;
			for(int x = 0; x < exist.length; x++){

				double[] temp = ArrayUtils.toPrimitive(exist[x]); 	
				DescriptiveStatistics desc = new DescriptiveStatistics(temp);
				double max = desc.getMax();
				double min = desc.getMin();
				if(all_max < max){

					all_max = max;
				}

				if(all_min > min){

					all_min = min;
				}

			}

			for(int x = 0; x < exist.length; x++){

				double[] temp = ArrayUtils.toPrimitive(exist[x]); 

				for(int y =0; y < temp.length; y++){
					norm_obj[x][y] = new Double((temp[y]-all_min)/(all_max-all_min));
				}
				System.out.println("profile_pos_norm_matrix:" + keys_list.get(q) + ":pos." + x + ":" + Arrays.toString(norm_obj[x]));
			}
			profile_pos_norm_matrix.put(keys_list.get(q),norm_obj);
		}//q
	}//method

	LinkedHashMap<String,Double[][]> profile_uniq_ligand;//Double[][]: position.uniq_ligand
	LinkedHashMap<String,Double[][]> profile_ligand_name;//Double[][]: position.ligand_name_indi

	public void map_profile_ligand_name() throws IOException{

	 	profile_ligand_name = new LinkedHashMap();
		profile_uniq_ligand = new LinkedHashMap();
	 	
	  	Pattern space_pattern = Pattern.compile("\\s+");
		Pattern tab_pattern = Pattern.compile("\\t+");
		Pattern comma_pattern = Pattern.compile(",");
		Pattern colon_pattern = Pattern.compile(":");
		
		String input="art_foundation_log_online/data/nbdb.profile_ligand.txt.3";
		
		List<String> lines  = FileUtils.readLines(new File(input));
		
		//Profile=ADEP
		//N	NAP_N7N	0	0	0	0	0	0	0	0	0	0	0	0	0	0	11	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
		//B	NAD_N6A	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0

		for(int q= 0; q < lines.size(); q++){
		
			String one_line = lines.get(q).trim();
			
			
			if(one_line.startsWith("Profile=")){
			
				String pro_name = one_line.replaceFirst("Profile=","").trim();
				System.out.println(pro_name);
				
				for(int w= q+1; w < lines.size(); w++){
		
					String next_line = lines.get(w).trim();

					if(!next_line.equals("")){
					
						if(next_line.startsWith("Profile=")){

							break;
						}

						
						String[] space_split = space_pattern.split(next_line);
						if(space_split.length == 1){
							System.out.println(next_line);
						}
						
						String uniq_temp = space_split[0];
						String ligand_temp = space_split[1];
						System.out.println("ligand name ind:" + ligand_temp);
						double[] freq = new double[space_split.length-2];

						for(int i = 2; i < space_split.length; i++){
							freq[i-2] = new Double(space_split[i]).doubleValue();

						}

						//ligand name start
						//if(ligand_name_list.contains(ligand_temp)){

							int ligand_name_ind = -1;

							for(int x = 0; x < ligand_name_indi.length; x++){
								if(ligand_temp.toUpperCase().contains(ligand_name_indi[x].toUpperCase())){
									ligand_name_ind = x;
									break;

								}

							}
							System.out.println("ligand name ind:" + ligand_name_ind);

							if(ligand_name_ind != -1){
								if(profile_ligand_name.containsKey(pro_name)){

									Double[][] exist = profile_ligand_name.get(pro_name);

									if(freq.length == exist.length){
										for(int x = 0; x < exist.length; x++){

											exist[x][ligand_name_ind] = exist[x][ligand_name_ind]  +new Double(freq[x]);
													
											
										}
										profile_ligand_name.put(pro_name,exist);
									}else{

										System.out.println("error length not same 499:" + ligand_temp + ":" + next_line + ":" + exist.length);
									}

								}else{
						
									Double[][] exist = new Double[freq.length][ligand_name_indi.length];
									for(int x = 0; x < freq.length; x++){

										for(int y = 0; y < ligand_name_indi.length; y++){
											
											if(freq[x] != 0 && y == ligand_name_ind){

												exist[x][y] = new Double(freq[x]);
											}else{

												exist[x][y] = new Double(0.0);
											}
										}
										
									}
									profile_ligand_name.put(pro_name,exist);
								}//else
							}else{
								System.out.println("ligand name error:" + next_line);

							}
						//}

						///uniq ligand start
						if(uniq_ligand_list.contains(uniq_temp)){

							int uniq_ligand_ind = -1;

							for(int x = 0; x < uniq_ligand.length; x++){
								if(uniq_temp.toUpperCase().contains(uniq_ligand[x].toUpperCase())){
									uniq_ligand_ind = x;
									break;

								}

							}

							if(profile_uniq_ligand.containsKey(pro_name)){

								Double[][] exist = profile_uniq_ligand.get(pro_name);

								if(freq.length == exist.length){
									for(int x = 0; x < exist.length; x++){

										exist[x][uniq_ligand_ind] = exist[x][uniq_ligand_ind]  +new Double(freq[x]);
												
										
									}
									profile_uniq_ligand.put(pro_name,exist);
								}else{

									System.out.println("error length not same 499:" + uniq_temp + ":" + next_line + ":" + exist.length);
								}

							}else{
					
								Double[][] exist = new Double[freq.length][uniq_ligand.length];
								for(int x = 0; x < freq.length; x++){

									for(int y = 0; y < uniq_ligand.length; y++){
										
										if(freq[x] != 0 && y == uniq_ligand_ind){

											exist[x][y] = new Double(freq[x]);
										}else{

											exist[x][y] = new Double(0.0);
										}
									}
									
								}
								profile_uniq_ligand.put(pro_name,exist);
							}//else
						}
					}//not ""
				}//w
					
			}//startswith		
					
		}//q

		List<String> keys_list = new LinkedList<String>(profile_uniq_ligand.keySet());	

		for(int q = 0; q < keys_list.size();	q++){

			Double[][] exist = profile_uniq_ligand.get(keys_list.get(q));
			
			for(int x = 0; x < exist.length; x++){

				 	
				System.out.println("profile_uniq_ligand:" + keys_list.get(q) + ":pos." + x + ":" + Arrays.toString(exist[x]));
			}
		}//q

		keys_list = new LinkedList<String>(profile_ligand_name.keySet());	

		for(int q = 0; q < keys_list.size();	q++){

			Double[][] exist = profile_ligand_name.get(keys_list.get(q));
			
			for(int x = 0; x < exist.length; x++){

				 	
				System.out.println("profile_ligand_name:" + keys_list.get(q) + ":pos." + x + ":" + Arrays.toString(exist[x]));
			}
		}//q

	//System.out.println("profile_ligand_name.size():"+profile_ligand_name.size());
	}//method
}//NBDB class


}//class
