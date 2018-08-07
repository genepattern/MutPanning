/************************************************************           
 * MutPanning - Step 12										*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: This is script is the first step to filter out	*
 * false positive significant genes. In this step the script*
 * searches for "hotspot" positions in significantly mutant	*
 * genes and writes the context around these hotspots into	*
 * a fastq file for blat query . The purpose of this script	*
 * is to determine read misalignments in highly repetitive 	*
 * and low complex regions as a potential source for 		*
 * deviation from the background mutation 					*
 * distribution pattern.									*
 *************************************************************/


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;



public class Filter_Step1 {
	
	static String[] chr={"1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","22","X","Y"};	
	static ArrayList<Gene>[] genes=new ArrayList[chr.length];
	
	static String[] entities=new String[0];
	static String[] files_sign=new String[2];
	static String file_genes="";
	static String[] file_out=new String[2];
	static String file_count="";
	static String[] index_header_samples={"ID","Sample","Cohort"};
	
	
	static String[] exclude={
			"PPIAL4A",//too long
			"PPIAL4B",
			"PPIAL4C",
			"FAM75A5",
			"FAM75A7",
			"ANKRD20A3",
			"FOXD4L2",
			"FOXD4L4",
			"FAM21B",
			"ASMT",//twice
			"ASMTL",
			"CSF2RA",
			"DHRSX",
			"IL3RA",
			"IL3R",
			"IL9R",
			"SPRY3",
			"ZBED1"};
	
	/*
	 * argument0: root file
	 * argument1: file samples
	* argument2: path to the Hg19 folder
	 */
	
	public static void main(String[] args){
		
		files_sign=new String[]{args[0]+"SignificanceRaw/Significance",args[0]+"SignificanceRaw/SignificanceUniform"};
		file_genes=args[2]+”Exons_Hg19.txt";
		file_out=new String[]{args[0]+"PostSignFilter/Queries/Query",args[0]+"PostSignFilter/Queries/QueryUniform"};
		file_count=args[0]+"EntityCounts/EntityCounts";
		
		if(!new File(args[0]+"PostSignFilter/Queries/").exists()){
			new File(args[0]+"PostSignFilter/Queries/").mkdirs();
		}
		
	
		
		//Determine all entity names, as clustering is performed separately for each Step this is needed to coordinate the order
		try{
			FileInputStream in=new FileInputStream(args[1]);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			ArrayList<String> aa =new ArrayList<String>();
			while((s=input.readLine())!=null){
				String e=s.split("	")[index_header[2]];
				if(!contains(e,aa)){
					aa.add(e);
				}
			}
			input.close();
			Collections.sort(aa);
			aa.add("PanCancer");
			entities=new String[aa.size()];
			for (int i=0;i<aa.size();i++){
				entities[i]=aa.get(i);
			}
		}
		catch(Exception e){
			System.out.println(e);
		}		
		
		
		try{
			
			for (int i=0;i<genes.length;i++){
				genes[i]=new ArrayList<Gene>();
			}
			String s="";
			
			//read the coordinates of all genes
			FileInputStream in=new FileInputStream(file_genes);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				if(index(t[0],exclude)!=-1){
					continue;
				}
				int ii=index_gene(t[0],genes[Integer.parseInt(t[1])-1]);
				if(ii==-1){
					genes[Integer.parseInt(t[1])-1].add(new Gene(t[0],Integer.parseInt(t[2]),Integer.parseInt(t[3])));
				}
				else{
					genes[Integer.parseInt(t[1])-1].get(ii).coord.add(new int[]{Integer.parseInt(t[2]),Integer.parseInt(t[3])});
				}
			}
			input.close();
			for (int i=0;i<genes.length;i++){
				for (int j=0;j<genes[i].size();j++){
					genes[i].get(j).start=min(genes[i].get(j).coord);
					genes[i].get(j).end=max(genes[i].get(j).coord);
				}
			}
			
			Comparator<Gene> comp=(Gene a, Gene b)->{
				return new Integer(a.start).compareTo(new Integer(b.start));
			};
			
			for (int i=0;i<genes.length;i++){
				Collections.sort(genes[i],comp);
			}
			
			for (int l=0;l<2;l++){
				for (int k=0;k<entities.length;k++){
					if(!new File(files_sign[l]+entities[k]+".txt").exists()){
						continue;
					}
					
					
					//extract all significantly mutant genes. these are the regions in which the
					//script searches for hotspot positions
					in=new FileInputStream(files_sign[l]+entities[k]+".txt");
					inn=new DataInputStream(in);
					input= new BufferedReader(new InputStreamReader(inn));
					input.readLine();
				
					ArrayList<String> sign_genes=new ArrayList<String>();
					while((s=input.readLine())!=null){
						String[] t=s.split("	");
						if(Double.parseDouble(t[15])<=0.25){
							sign_genes.add(t[0]);
						}
					}
					input.close();
					
					//handle the search for hotspot postions separately for each chr
					Subthread[] threads=new Subthread[chr.length];
					for (int i=0;i<threads.length;i++){
						threads[i]=new Subthread();
						threads[i].sign_genes=sign_genes;
						threads[i].entity=entities[k];
						threads[i].c=i;
						threads[i].start();

					}

					//wait until the computations are done for each chr
					boolean all_done=false;
					do{
						Thread.sleep(3000);
						all_done=true;
						for (int i=0;i<threads.length;i++){
							if(!threads[i].done){
								all_done=false;
								break;
							}
						}
					}while(!all_done);
					
					
					//output of the query in the format of an fastq file
					FileWriter out=new FileWriter(file_out[l]+entities[k]+".fa");
					BufferedWriter output= new BufferedWriter(out);
					
					for (int i=0;i<threads.length;i++){
						for (int j=0;j<threads[i].names.size();j++){
							output.write(threads[i].names.get(j));
							output.newLine();
							output.write(threads[i].seq.get(j));
							output.newLine();
						}
					}
					output.close();
				}
			}
			
		}
		catch(Exception e){
			System.out.println(e);
		}
	}
	
	public static int[] index_header(String[] header, String[] ideal_header){
		int[] indices=new int[ideal_header.length];
		for (int i=0;i<ideal_header.length;i++){
			int index=-1;
			for (int j=0;j<header.length;j++){
				if(header[j].equals(ideal_header[i])){
					index=j;
					break;
				}
			}
			indices[i]=index;
		}
		return indices;
	}
	
	public static boolean contains(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return true;
			}
		}
		return false;
	}
	
	private static class Gene{
		String name="";
		ArrayList<int[]> coord=new ArrayList<int[]>();
		int start=-1;
		int end=-1;
		
		
		public Gene(String name, int start, int end){
			this.name=name;
			coord.add(new int[]{start,end});
		}
		public boolean contains (int pos){
			for (int i=0;i<coord.size();i++){
				if(coord.get(i)[0]<=pos&&pos<=coord.get(i)[1]){
					return true;
				}
			}
			return false;
		}
	}
	
	public static int min(ArrayList<int[]> a){
		int min=1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[0]<min){
				min=a.get(i)[0];
			}
		}
		return min;
	}
	
	public static int max(ArrayList<int[]> a){
		int max=-1000000000;
		for (int i=0;i<a.size();i++){
			if(a.get(i)[1]>max){
				max=a.get(i)[1];
			}
		}
		return max;
	}
	
	public static int index(String s, ArrayList<String> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).equals(s)){
				return i;
			}
		}
		return -1;
	}
	public static int index_gene(String s, ArrayList<Gene> t){
		for (int i=0;i<t.size();i++){
			if(t.get(i).name.equals(s)){
				return i;
			}
		}
		return -1;
	}
	public static int index(String s, String[] t){
		for (int i=0;i<t.length;i++){
			if(t[i].equals(s)){
				return i;
			}
		}
		return -1;
	}
	
	private static class Subthread extends Thread{
		volatile boolean done=false;
		int c=-1;
		ArrayList<String> names=new ArrayList<String>();
		ArrayList<String> seq=new ArrayList<String>();
		ArrayList<String> sign_genes=new ArrayList<String>();
		String entity="";
		public void run(){
			System.out.println("start "+chr[c]);
			
			try{
				FileInputStream in=new FileInputStream(file_count+entity+"_Chr"+chr[c]+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				
				ArrayList<Integer> position=new ArrayList<Integer>();
				ArrayList<double[]> lambda=new ArrayList<double[]>();
				ArrayList<String> nucl=new ArrayList<String>();
				ArrayList<Double>coverage=new ArrayList<Double>();
				ArrayList<int[]> label=new ArrayList<int[]>();
				ArrayList<int[]> count=new ArrayList<int[]>();
				
				
				//go through the reference squences of the genes and search for hotspot positions (at least 4 mutations)
				//if the gene appears in the list of sign genes the sequence context aorund each
				//hotspot position is saved as a fastq format. to have the context around a hotspot postions the reference
				//sequence of each gene is saved in a queue.
				for (int nn=0;nn<genes[c].size();nn++){
					int ii=0;
					if(position.size()>0){
						while(ii<position.size()&&position.get(ii)<genes[c].get(nn).start){
							ii++;
						}
					}
					if(ii>0){
						for (int i=ii-1;i>=0;i--){
							nucl.remove(i);
							lambda.remove(i);
							position.remove(i);
							coverage.remove(i);
							count.remove(i);
							label.remove(i);
						}
					}
					
					String s="";
					while ((s=input.readLine())!=null){
						String[] t=s.split("	");
						position.add(Integer.parseInt(t[0]));
						nucl.add(t[1]);
						coverage.add(Double.parseDouble(t[2]));
						String[] tt=t[3].split(";");
						label.add(new int[]{Integer.parseInt(tt[0]),Integer.parseInt(tt[1]),Integer.parseInt(tt[2])});
						tt=t[4].split(";");
						lambda.add(new double[]{Double.parseDouble(tt[0]),Double.parseDouble(tt[1]),Double.parseDouble(tt[2])});
						tt=t[5].split(";");
						count.add(new int[]{Integer.parseInt(tt[0]),Integer.parseInt(tt[1]),Integer.parseInt(tt[2])});
					
						if(Integer.parseInt(t[0])>=genes[c].get(nn).end){
							break;
						}
					}

					if(index(genes[c].get(nn).name,sign_genes)!=-1){
						int counter=1;
						for (int i=0;i<position.size();i++){
							if(!genes[c].get(nn).contains(position.get(i))){
								continue;
							}
							for (int j=0;j<3;j++){
								if(count.get(i)[j]<4){
									continue;
								}
								String ref_seq="";
								String alt_seq="";
								int max_index=-100;
								int min_index=-100;
								for (int k=-25;k<=25;k++){
									if(k+i>=0&&k+i<position.size()&&position.get(k+i)==position.get(i)+k){
										max_index=k;
										if(min_index==-100){
											min_index=k;
										}
										if(k!=0){
											ref_seq=ref_seq+nucl.get(i+k);
											alt_seq=alt_seq+nucl.get(i+k);
										}
										else{
											ref_seq=ref_seq+nucl.get(i+k);
											alt_seq=alt_seq+mut(nucl.get(i+k),j);
										}
									}
								}
								//names.add(">"+genes[c].get(nn).name+"_ref_"+counter+";"+min_index+","+max_index+";"+count.get(i)[j]+";"+chr[c]+","+position.get(i));
								//seq.add(ref_seq);
								names.add(">"+genes[c].get(nn).name+"_alt_"+counter+";"+min_index+","+max_index+";"+count.get(i)[j]+";"+chr[c]+","+position.get(i)+","+j);
								seq.add(alt_seq);
								counter++;
							}
							
						}
					}
					
				}
				input.close();
			}
			catch(Exception e){
				System.out.println(e);
			}
			System.out.println("done "+chr[c]);
			done=true;
		}
	}
	
	
	public static String mut(String nucl, int c){
		if(nucl.equals("A")){
			return new String[]{"T","C","G"}[c];
		}
		else if(nucl.equals("C")){
			return new String[]{"A","G","T"}[c];
		}
		else if(nucl.equals("G")){
			return new String[]{"T","C","A"}[c];
		}
		else if(nucl.equals("T")){
			return new String[]{"A","G","C"}[c];
		}
		else{
			return "N";
		}
	}
}
