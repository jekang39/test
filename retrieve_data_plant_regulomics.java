
import java.io.*;
import java.util.*;
import java.util.regex.*;
import org.apache.commons.io.IOUtils;
import org.apache.commons.io.FileUtils;

public class  retrieve_data_plant_regulomics {


public static void main(String[] args) throws IOException,InterruptedException{

	retrieve_data_plant_regulomics aa = new retrieve_data_plant_regulomics();
	aa.map_pfam_clans_only();
	aa.map_pfam_arabi();
	aa.read_go_po_files();
	aa.get_plant_regulomics_ppi_pfam_assign();
	aa.get_plant_regulomics_motif_result();
	aa.read_iregnet_result();
	
}

LinkedHashMap<String, String> pfam_clan_map_only;


void map_pfam_clans_only() throws IOException{

	//map pfam to clan. clan is sort of superfamily of pfam
	System.out.println("map_pfam_clans_only():" );

	pfam_clan_map_only = new LinkedHashMap();
 	
 	
  	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	
	String input="data/Pfam-A.clans.tsv";
	
	List<String> lines  = FileUtils.readLines(new File(input));
	
	for(int q= 0; q < lines.size(); q++){
	
		String one_line = lines.get(q).trim().toUpperCase();
		
		if(!one_line.equals("")){
		
			String[] one_line_tab = tab_pattern.split(one_line);
			
			if(one_line_tab.length != 5){
			
				//System.out.println("tab error:"+one_line_tab.length + ":" +q + ":"+one_line);
			}
			String clan ="";
			if(one_line_tab[1].trim().equals("NA")){

				clan = "NA";
			
			}else{
			
				clan = one_line_tab[1].toUpperCase();
				
			
			}
			String pfamid = one_line_tab[0].toUpperCase();
			
			if(!pfam_clan_map_only.containsKey(pfamid)){
				
				//pfam_clan_map_only.put(pfamid, one_line_tab[1]+":"+one_line_tab[2]+":"+one_line_tab[3]+":"+one_line_tab[4]);
				pfam_clan_map_only.put(pfamid, clan);	
				
						
			}else{
				//System.out.println(q + ":double key error:" + pfamid+":"+one_line);
			}
			
			
			
		}//if
			

}//q
System.out.println("pfam_clan_map_only:" + pfam_clan_map_only.size());
}

LinkedHashMap<String,LinkedList<String>> arabi_to_pfam;


void map_pfam_arabi() throws IOException{

	//per arabidopsis gene, atxg, map all pfams assigned to the gene

	Pattern colon_pattern = Pattern.compile(":");
	Pattern tab_pattern = Pattern.compile("\\t");
	Pattern comma_pattern = Pattern.compile(",");

	arabi_to_pfam = new LinkedHashMap();

	
	String path = "data/map.atxg.pfam";

    List<String> lines  = FileUtils.readLines(new File(path));
	
	for(int q= 0; q < lines.size(); q++){
	
	
		String one_line = lines.get(q).trim();
		String[] split_str = colon_pattern.split(one_line);
		
		String id = split_str[0].trim();
		
		String pfams = split_str[1].trim();
		String[] split_pfam = comma_pattern.split(pfams);
		LinkedList<String> exist = new LinkedList();
		for(int i = 0; i < split_pfam.length; i++){

			exist.add(split_pfam[i].trim());
		}
		
		
		arabi_to_pfam.put(id,exist);
		
			
 	}//q
 	
 	
 	System.out.println("arabi_to_pfam:" +arabi_to_pfam.size());

	List<String> keys_list = new LinkedList<String>(arabi_to_pfam.keySet());
	
				
	for(int k = 0; k < keys_list.size(); k++){
				
		String id = keys_list.get(k).trim();

		LinkedList<String> exist = arabi_to_pfam.get(id);
		System.out.println("arabi_to_pfam:" + id + ":"  + Arrays.toString(exist.toArray()));
		
	}

}

LinkedHashMap<String,LinkedList<String>>[] map_go_po_pfam;
String[] input = {"pfam", "signal","slim","asso","ana","tempo","feature"};


 void read_go_po_files() throws IOException{
	
//file contains genes containingWD40 and those interact with genes containing wd40- results from iregnet. for these genes, get pfam assignment, go-po annotation

	String[] input_files = {"data/plant_regulomics/out_map_arabi_pfam_arabidopsis.pfam.assign.wd40.iregnet","data/plant_regulomics/out_map_arabi_pfam_arabidopsis.signal.wd40.iregnet","data/plant_regulomics/out_map_arabi_pfam_arabidopsis.slim.wd40.iregnet","data/plant_regulomics/out_map_arabi_pfam_arabidopsis.asso.wd40.iregnet","data/plant_regulomics/out_map_arabi_pfam_arabidopsis.ana.wd40.iregnet","data/plant_regulomics/out_map_arabi_pfam_arabidopsis.tempo.wd40.iregnet","data/plant_regulomics/out_map_arabi_pfam_arabidopsis.feature.wd40.iregnet"};

	String[] remove_str = {"NANANA","arabi_to_signal:","arabi_to_slim:","arabi_to_asso:","arabi_to_ana:","arabi_to_tempo:","arabi_to_feature:"};
 
 	System.out.println("read_go_po_files():" );
 
 	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");
	map_go_po_pfam = new LinkedHashMap[input_files.length];

	for(int j =0; j < input_files.length; j++){

		map_go_po_pfam[j] = new LinkedHashMap();

	 	List<String> lines  = FileUtils.readLines(new File(input_files[j]));
						
		for(int k = 0; k < lines.size(); k++){

			String one_line = lines.get(k);
			one_line = one_line.replaceFirst(remove_str[j], "");
			String id = one_line.substring(0,one_line.indexOf(":"));
			String inter_list = one_line.substring(one_line.indexOf(":")+1, one_line.length());
			inter_list = inter_list.replaceFirst("\\[","");
			inter_list = inter_list.replaceFirst("\\]","");
			String[] comma_split = comma_pattern.split(inter_list);
			LinkedList<String> list = new LinkedList();
			for(int i =0; i < comma_split.length; i++){
				String temp = comma_split[i].trim();
				list.add(temp);
				
			}


			//System.out.println("check:" + input[j] + ":" +id + ":" + Arrays.toString(list.toArray()));
			
			map_go_po_pfam[j].put(id,list);

		}
	}//j
}//method

void get_plant_regulomics_ppi_pfam_assign() throws IOException{
//for gene of interest, we get ppi. with all proteins that interact with the gene of interest, we retrieve pfams and make pfam list. loop over pfam list --> convert to clan list where clan is sort of superfamily of pfam. this fuction prints gene of interest and clan pool created based on ppi
 
 	System.out.println("get_plant_regulomics_ppi_pfam_assign():" );
	
	LinkedList<String> all_clan = new LinkedList();

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	String[] input_files = {"data/plant_regulomics/agl15.AT5G13790/AT5G13790_PPI.csv","data/plant_regulomics/lec1.AT1G21970/AT1G21970_PPI.csv","data/plant_regulomics/abi3.AT3G24650/AT3G24650_PPI.csv","data/plant_regulomics/fus3.AT3G26790/AT3G26790_PPI.csv","data/plant_regulomics/HSI2.AT2G30470/AT2G30470_PPI.csv","data/plant_regulomics/ARF1.at1g23490/AT1G23490_PPI.csv","data/plant_regulomics/ARF1.at1g59750/AT1G59750_PPI.csv","data/plant_regulomics/ARF5.at1g19850/AT1G19850_PPI.csv","data/plant_regulomics/MEA.AT1G02580/AT1G02580_PPI.csv","data/plant_regulomics/CLF.AT2G23380/AT2G23380_PPI.csv","data/plant_regulomics/SWN.AT4G02020/AT4G02020_PPI.csv"};

	String[] genes = {"AT5G13790","AT1G21970","AT3G24650","AT3G26790","AT2G30470","AT1G23490","AT1G59750","AT1G19850","AT1G02580","AT2G23380","AT4G02020"};

 	LinkedList<String>[] pfams_ppi = new LinkedList[input_files.length];
	LinkedList<String>[] clans_ppi = new LinkedList[input_files.length];

	for(int j =0; j < input_files.length; j++){

		pfams_ppi[j] = new LinkedList();
		clans_ppi[j] = new LinkedList();
		
		List<String> lines  = FileUtils.readLines(new File(input_files[j]));
						
		for(int k = 1; k < lines.size(); k++){

			String one_line = lines.get(k);
			
			String[] space_split = space_pattern.split(one_line);
			String id = space_split[0];
			LinkedList<String> pfams = arabi_to_pfam.get(id);
			if(pfams != null){
				pfams_ppi[j].addAll(pfams);
				//System.out.println("check:" + genes[j] + ":" +id + ":" + Arrays.toString(pfams.toArray()));
			}
			

		}

		Set<String> pfams_h_inv = new LinkedHashSet<>();
		pfams_h_inv.addAll((List)pfams_ppi[j]);
										
		// Clear the list
		((List)pfams_ppi[j]).clear();
																  
		// add the elements of set
		// with no duplicates to the list
		((List)pfams_ppi[j]).addAll(pfams_h_inv);
		System.out.println("all pfam:" + genes[j] + Arrays.toString(pfams_ppi[j].toArray()));

		for(int i = 0; i < pfams_ppi[j].size(); i++){
		
			String clan = pfam_clan_map_only.get(pfams_ppi[j].get(i));

			if(!clans_ppi[j].contains(clan)){

				clans_ppi[j].add(clan);
			}

			if(!all_clan.contains(clan)){

				all_clan.add(clan);
			}

		}
		
	}//j

	System.out.println("all_clan:" + Arrays.toString(all_clan.toArray()));

	for(int j =0; j < input_files.length; j++){

		String write = genes[j] + ":clan:";

		for(int i = 0; i < all_clan.size(); i++){
		
			String clan = all_clan.get(i);

			if(clans_ppi[j].contains(clan)){

				write = write + clan + ":";
			}else{

				write = write +  "NA:";
			}

		}

		System.out.println(write);
	}

}//method	

void get_plant_regulomics_motif_result() throws IOException{
 //without GRN construction, it is confusing. here, I combined all information (gene of interest, its binding factor, those need to be compared with gene of interest, and those need to be compared with binding factor etc.) we only retrieve information that may interact. we put all genes in the same array. It is beyond the scope of this article. in this article, we only try to show typical scenario of PNI/PPI where ART FOUNDATION-LOG can be used and peach homologues- how linked to peach phenotypes. In biomarker discovery with gwas data which contain a large number of genetic variants, we need to construct GRN, and integrate with EWAS, TWAS, GWAS in applying ART FOUNDATION-LOG. it remains for further studies

 	System.out.println("get_plant_regulomics_motif_result():" );

	LinkedList<String> all_promoter_class = new LinkedList();
	LinkedList<String> all_gene_body_class = new LinkedList();

	String[] header = {"Name","Chr","Start","End","Sequence","Motif_Score","JasparID","Specie","Class"};
	String[] genes = {"AT5G13790","AT1G21970","AT3G24650","AT3G26790","AT2G30470","AT1G23490","AT1G59750","AT1G19850","AT1G02580","AT2G23380","AT4G02020"};
	String[] strand = {"-", "-","+","-","-","+","+","+","+","+","+"};
	int[] tss_five_prime_utr_start = {4450843,7729617,8997370,9855989,12985043,8336683,21979384,6886879,544783,9955553,886600};
	

	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	String[] input_files = {"data/plant_regulomics/agl15.AT5G13790/AT5G13790_motif_result.csv","data/plant_regulomics/lec1.AT1G21970/AT1G21970_motif_result.csv","data/plant_regulomics/abi3.AT3G24650/AT3G24650_motif_result.csv","data/plant_regulomics/fus3.AT3G26790/AT3G26790_motif_result.csv","data/plant_regulomics/HSI2.AT2G30470/AT2G30470_motif_result.csv","data/plant_regulomics/ARF1.at1g23490/AT1G23490_motif_result.csv","data/plant_regulomics/ARF1.at1g59750/AT1G59750_motif_result.csv","data/plant_regulomics/ARF5.at1g19850/AT1G19850_motif_result.csv","data/plant_regulomics/MEA.AT1G02580/AT1G02580_motif_result.csv","data/plant_regulomics/CLF.AT2G23380/AT2G23380_motif_result.csv","data/plant_regulomics/SWN.AT4G02020/AT4G02020_motif_result.csv"};

 	LinkedHashMap<String,LinkedList<String>>[] promoter_family_gene = new LinkedHashMap[input_files.length];
	LinkedHashMap<String,LinkedList<Integer[]>>[] promoter_family_index = new LinkedHashMap[input_files.length];
	LinkedHashMap<String,LinkedList<String>>[] promoter_family_motif = new LinkedHashMap[input_files.length];
	LinkedHashMap<String,LinkedList<String>>[] gene_body_family_gene = new LinkedHashMap[input_files.length];
	LinkedHashMap<String,LinkedList<Integer[]>>[] gene_body_family_index = new LinkedHashMap[input_files.length];
	LinkedHashMap<String,LinkedList<String>>[] gene_body_family_motif = new LinkedHashMap[input_files.length];

	for(int j =0; j < input_files.length; j++){

		promoter_family_gene[j] = new LinkedHashMap();
		promoter_family_motif[j] = new LinkedHashMap();
		promoter_family_index[j] = new LinkedHashMap();
		gene_body_family_gene[j] = new LinkedHashMap();
		gene_body_family_index[j] = new LinkedHashMap();
		gene_body_family_motif[j] = new LinkedHashMap();
		
		List<String> lines  = FileUtils.readLines(new File(input_files[j]));
						
		for(int k = 1; k < lines.size(); k++){

			String one_line = lines.get(k);
			
			String[] tab_split = tab_pattern.split(one_line);
			String id = tab_split[0];
			String motif = tab_split[6];
			String family = tab_split[8];
			Integer start = new Integer(-1);
			Integer end = new Integer(-1);
			boolean promoter_flag = false;

			if(strand[j].equals("+")){
		
				start = new Integer(tab_split[2]);
				if(start.intValue() < tss_five_prime_utr_start[j]){

					promoter_flag = true;
				}
				
				end = new Integer(tab_split[3]);

			}else{

				int start_int = new Integer(tab_split[3]).intValue();
				start_int = Math.abs(start_int);
				start = new Integer(start_int);
				if(start_int > tss_five_prime_utr_start[j]){

					promoter_flag = true;
				}
				int end_int = new Integer(tab_split[2]).intValue();
				end_int = Math.abs(end_int);
				end = new Integer(end_int);
			}

			Integer[] temp = new Integer[2];
			temp[0] = start;
			temp[1] = end;

			if(promoter_flag){

				if(promoter_family_gene[j].containsKey(family)){
			
					LinkedList<String> exist = promoter_family_gene[j].get(family);
					exist.add(id);
					
					promoter_family_gene[j].put(family,exist);
					
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(id);
					promoter_family_gene[j].put(family,exist);
				}

				if(promoter_family_index[j].containsKey(family)){
			
					LinkedList<Integer[]> exist = promoter_family_index[j].get(family);
					
					exist.add(temp);
					
					promoter_family_index[j].put(family,exist);
					
				}else{
					LinkedList<Integer[]> exist = new LinkedList();

					exist.add(temp);
					promoter_family_index[j].put(family,exist);
				}

				if(promoter_family_motif[j].containsKey(family)){
			
					LinkedList<String> exist = promoter_family_motif[j].get(family);
					exist.add(motif);
					
					promoter_family_motif[j].put(family,exist);
					
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(motif);
					promoter_family_motif[j].put(family,exist);
				}


			}else{

				if(gene_body_family_gene[j].containsKey(family)){
			
					LinkedList<String> exist = gene_body_family_gene[j].get(family);
					exist.add(id);
					
					gene_body_family_gene[j].put(family,exist);
					
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(id);
					gene_body_family_gene[j].put(family,exist);
				}

				if(gene_body_family_index[j].containsKey(family)){
			
					LinkedList<Integer[]> exist = gene_body_family_index[j].get(family);
					
					exist.add(temp);					
					gene_body_family_index[j].put(family,exist);
					
				}else{
					LinkedList<Integer[]> exist = new LinkedList();

					exist.add(temp);
					gene_body_family_index[j].put(family,exist);
				}

				if(gene_body_family_motif[j].containsKey(family)){
			
					LinkedList<String> exist = gene_body_family_motif[j].get(family);
					exist.add(motif);
					
					gene_body_family_motif[j].put(family,exist);
					
				}else{
					LinkedList<String> exist = new LinkedList();
					exist.add(motif);
					gene_body_family_motif[j].put(family,exist);
				}

			}

		}//k
	}//j

	List<String>[] promoter_family_list = new List[genes.length];
	List<String>[] gene_body_family_list = new List[genes.length];

	for(int j =0; j < genes.length; j++){

		promoter_family_list[j] = new LinkedList<String>(promoter_family_gene[j].keySet());
		System.out.println(genes[j] + ":promoter - all motif Classes:" + Arrays.toString(promoter_family_list[j].toArray()));

		//make all promoter motif class

		all_promoter_class.addAll(promoter_family_list[j]);

		for(int q = 0; q < promoter_family_list[j].size(); q++){
		
			String family = promoter_family_list[j].get(q);
			LinkedList<String> ids = promoter_family_gene[j].get(family);
			System.out.println(genes[j] + ":promoter motif Class:" + family + ":genes:" +  Arrays.toString(ids.toArray()));
			
			
		}

		gene_body_family_list[j]  = new LinkedList<String>(gene_body_family_gene[j].keySet());
		System.out.println(genes[j] + ":gene_body - all motif Classes:" + Arrays.toString(gene_body_family_list[j].toArray()));

		//make all gene_body motif class

		all_gene_body_class.addAll(gene_body_family_list[j]);

		for(int q = 0; q < gene_body_family_list[j].size(); q++){
		
			String family = gene_body_family_list[j].get(q);
			LinkedList<String> ids = gene_body_family_gene[j].get(family);
			System.out.println(genes[j] + ":gene_body motif Class:" + family + ":genes:" +  Arrays.toString(ids.toArray()));
			LinkedList<Integer[]> pos = gene_body_family_index[j].get(family);
			
		}
		
	}//j

	//get all promoter classes

	Set<String> pfams_h_inv = new LinkedHashSet<>();
	pfams_h_inv.addAll((List)all_promoter_class);
											
	// Clear the list
	((List)all_promoter_class).clear();
																	  
	// add the elements of set
	// with no duplicates to the list
	((List)all_promoter_class).addAll(pfams_h_inv);

	System.out.println("all_promoter_class:" + Arrays.toString(all_promoter_class.toArray()));

	for(int j =0; j < input_files.length; j++){

		String write = genes[j] + ":promoter_family:";

		for(int i = 0; i < all_promoter_class.size(); i++){
		
			String family = all_promoter_class.get(i);

			if(promoter_family_list[j].contains(family)){

				write = write + family + ":";
			}else{

				write = write +  "NA:";
			}

		}

		System.out.println(write);
	}

	//get all gene body classes

	pfams_h_inv = new LinkedHashSet<>();
	pfams_h_inv.addAll((List)all_gene_body_class);
											
	// Clear the list
	((List)all_gene_body_class).clear();
																	  
	// add the elements of set
	// with no duplicates to the list
	((List)all_gene_body_class).addAll(pfams_h_inv);
	System.out.println("all_gene_body_class:" + Arrays.toString(all_gene_body_class.toArray()));

	for(int j =0; j < input_files.length; j++){
		
		String write = genes[j] + ":gene_body family:";

		for(int i = 0; i < all_gene_body_class.size(); i++){
		
			String family = all_gene_body_class.get(i);

			if(gene_body_family_list[j].contains(family)){

				write = write + family + ":";
			}else{

				write = write +  "NA:";
			}

		}

		System.out.println(write);
	}


}//method	

LinkedHashMap<String,LinkedList<String>> map_iregnet_target;
LinkedHashMap<String,LinkedList<String>> map_iregnet_ppi;
 
 void read_iregnet_result() throws IOException{
 
 	System.out.println("read_iregnet_result():" );
 
 	
	Pattern space_pattern = Pattern.compile("\\s+");
	Pattern tab_pattern = Pattern.compile("\\t+");
	Pattern comma_pattern = Pattern.compile(",");
	Pattern colon_pattern = Pattern.compile(":");

	map_iregnet_target= new LinkedHashMap();
	map_iregnet_ppi= new LinkedHashMap();
	
 
 	String path = "data/plant_regulomics/mqCUz4ACBVBGJjckdfi9O3AX6HeQI4Xk.multi_gene.csv";
	List<String> lines  = FileUtils.readLines(new File(path));
					
	for(int k = 0; k < lines.size(); k++){

		/*
		Name	AGI ID	TF/Histone	P-value	Corr. p-value	GSE	Genotype	Tissue	Growth cond.	IP	COEX	COEX (Tissue)	COEX (Abiotic)	COEX (Hormone)	COEX (Biotic)	COEX (Light)	Targets	PPI	INDEX
607	APRR5, PRR5	AT5G24470	Protein	1.37e-44	2.2413200000000002e-41	GSE36361	pPRR5::FLAG-PRR5-GFP prr5	whole plant	10 hours after light on	ChIP	0.5903389830508475	0.5462288135593221	1.446137288135593	0.9240305084745762	0.9489474576271187	1.179484745762712	AT1G01520, AT1G01580, AT1G06000, AT1G10370, AT1G10740, AT1G10960, AT1G15980, AT1G19150, AT1G19490, AT1G24147, AT1G30530, AT1G32080, AT1G35560, AT1G48100, AT1G49630, AT1G55910, AT1G55960, AT1G56510, AT1G64780, AT1G78570, AT1G80340, AT2G04790, AT2G17780, AT2G27402, AT2G27420, AT2G30100, AT2G37040, AT2G39730, AT2G46830, AT2G46940, AT3G01060, AT3G01440, AT3G01550, AT3G02380, AT3G03770, AT3G10185, AT3G10420, AT3G13510, AT3G19450, AT3G22840, AT3G24520, AT3G30460, AT3G47070, AT3G48420, AT3G50560, AT3G51750, AT3G53830, AT3G55330, AT4G00490, AT4G01080, AT4G12830, AT4G14540, AT4G14690, AT4G27940, AT4G28660, AT4G30610, AT5G08640, AT5G10150, AT5G13170, AT5G13930, AT5G15948, AT5G17170, AT5G17300, AT5G17870, AT5G17890, AT5G24660, AT5G38410, AT5G38430, AT5G40950, AT5G49330	AT5G58350, AT3G04910, AT2G46790, AT2G46790, AT5G02810, AT5G02810, AT5G61380, AT5G61380, AT5G02810, AT5G02810, AT1G01060, AT2G46830, AT5G02810, AT5G02810, AT1G01060, AT2G46830, AT5G57360, AT5G57360, AT5G57360, AT5G57360, AT5G02810, AT5G02810, AT2G46790, AT2G46790, AT5G57360, AT5G57360, AT5G57360, AT2G18915, AT1G68050, AT5G57360, AT2G18915, AT1G68050, AT1G68050, AT5G61380, AT5G61380, AT5G61380, AT2G37000, AT1G69690, AT4G18390, AT1G15750, AT1G15750, AT1G15750, AT1G69690, AT1G60250, AT4G04885, AT4G18390, AT2G34450, AT4G39070, AT2G21320, AT5G24520, AT2G46670, AT2G18090, AT5G20240, AT3G47620, AT3G15030, AT4G34680, AT1G44830, AT4G38960	AT5G24470.Protein.GSE36361_0.ChIP
		*/

		String one_line = lines.get(k).toUpperCase();
		String[] tab_split = tab_pattern.split(one_line);
		String id = tab_split[2].trim();
		String target_list = tab_split[17].trim();
		
		String[] comma_split = comma_pattern.split(target_list);
		LinkedList<String> list = new LinkedList();
		for(int i =0; i < comma_split.length; i++){
			String temp = comma_split[i].trim();
			list.add(temp);
			
		}
		map_iregnet_target.put(id,list);
		target_list = tab_split[18].trim();
		
		comma_split = comma_pattern.split(target_list);
		list = new LinkedList();
		for(int i =0; i < comma_split.length; i++){
			String temp = comma_split[i].trim();
			list.add(temp);
			
		}
		map_iregnet_ppi.put(id,list);
	}
}//method

}//class
