/************************************************************           
 * MutPanning - Step 1									*		*
 * 															*   
 * Author:		Felix Dietlein								*   
 *															*   
 * Copyright:	(C) 2018 									*   
 *															*   
 * License:		Public Domain								*   
 *															*   
 * Summary: This scripts starts from a Maf file and 		*
 * 			annotates each position in the human			*
 * 			exome with the sample indices which carry		*
 * 			a mutation. These aligned files are needed		*
 * 			for the subsequent steps of the algorithm		*										*   
 *************************************************************/


import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Hashtable;

public class AlignHG19 {

	
	static String file_peptide="";
	static String file_coverage="";
	static String file_samples="";
	static String file_maf="";
	static String file_out="";
	static String[] chr={};
	static String[] index_header_samples={"ID","Sample","Cohort"};
	static String[] index_header_maf={"Hugo_Symbol","Chromosome","Start_Position","End_Position","Strand","Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Tumor_Sample_Barcode"};
	
	/*	argument 0: root file, where all the other files can be found
	 * argument 1: maf file
	 * argument 2: sample annotation file
	 * agrument 3: which chr should be analyzed
	* argument 4: path to the Hg19 folder
	 */
	
	public static void main(String[] args){
		
		file_peptide=args[4]+”ASAnnotationHg19/ASAnnotation_chr";
		file_coverage=args[4]+”CoverageExome_TCGA/Coverage_chr";
		file_samples=args[2];
		file_maf=args[1];
		file_out=args[0]+"AlignHg19/AlignHg19Chr";
		
		if(!new File(args[0]+"AlignHg19").exists()){
			new File(args[0]+"AlignHg19").mkdir();
		}
		
		if(Integer.parseInt(args[3])==0){
			chr=new String[]{"1","2","3"};
		}
		else if(Integer.parseInt(args[3])==1){
			chr=new String[]{"4","5","6","7","8","9"};
		}
		else if(Integer.parseInt(args[3])==2){
			chr=new String[]{"10","11","12","13","14","15","16","17"};
		}
		else if(Integer.parseInt(args[3])==3){
			chr=new String[]{"18","19","20","21","22","X","Y"};
		}
		else{
			System.exit(0);
		}
		
	
		try{
			Hashtable<String, Integer> sample_table=new Hashtable<String, Integer>();
			// read sample names and link them to indices
			FileInputStream in=new FileInputStream(file_samples);
			DataInputStream inn=new DataInputStream(in);
			BufferedReader input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header=index_header(input.readLine().split("	"),index_header_samples);
			String s="";
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				sample_table.put(t[index_header[1]],Integer.parseInt(t[index_header[0]]));
			}
			input.close();
			
			//read the sequence of all chr specified
			Subthread[] threads=new Subthread[chr.length];
			for (int i=0;i<threads.length;i++){
				threads[i]=new Subthread();
				threads[i].chr=chr[i];
				threads[i].start();
			}
			
			//wait until all the reading of all chr is done
			boolean all_done=false;
			do{
				Thread.sleep(10000);
				all_done=true;
				for (int i=0;i<threads.length;i++){
					if(!threads[i].done){
						all_done=false;
						break;
					}
				}
			}while(!all_done);
			
			//load all sequences in a big position table
			Hashtable<Integer,Integer>[] position_table=new Hashtable[chr.length];
			ArrayList<Position>[] position=new ArrayList[chr.length];
			for (int i=0;i<chr.length;i++){
				position_table[i]=threads[i].table;
				position[i]=threads[i].position;
			}
			
			
			//walk through the mutation annotation file. for each mutation link the index
			//of its sample to the genomic position of the reference genome
			in=new FileInputStream(file_maf);
			inn=new DataInputStream(in);
			input= new BufferedReader(new InputStreamReader(inn));
			int[] index_header_m=index_header(input.readLine().split("	"),index_header_maf);
			while((s=input.readLine())!=null){
				String[] t=s.split("	");
				
				int chr_index=index(t[index_header_m[1]],chr);
				if(chr_index==-1){
					continue;
				}
				int pos=Integer.parseInt(t[index_header_m[2]]);
				String ref=t[index_header_m[7]].toUpperCase();
				String tumor=t[index_header_m[9]].toUpperCase();
				if(tumor.equals(ref)||tumor.equals("")){
					tumor=t[index_header_m[8]].toUpperCase();
				}
				
				if(ref.length()!=1||tumor.length()!=1){
					continue;
				}
				if(!isNucleotide(ref)||!isNucleotide(tumor)){
					continue;
				}
				
				int sample_index=sample_table.get(t[index_header_m[10]]);
				int type=type(ref,tumor);
				if(type==-1){
					continue;
				}
				Integer ii=position_table[chr_index].get(pos);
				if(ii==null){
					continue;
				}
				
				if(position[chr_index].get(ii.intValue()).pos!=pos||!position[chr_index].get(ii.intValue()).nucl.equals(ref)){
					//System.out.println("Misalignment Hg19");
					//System.out.println(t[10]);
					continue;
				}
				position[chr_index].get(ii.intValue()).add(sample_index,type);
			}
			input.close();
			
			//now the position table contains all the positions together with the indices of the samples that contain mutations
			
			//output loop, generate a separate file for each chr
			for (int i=0;i<position.length;i++){
				FileWriter out=new FileWriter(file_out+chr[i]+".txt");
				BufferedWriter output= new BufferedWriter(out);
				
				for (int j=0;j<position[i].size();j++){
					output.write(position[i].get(j).pos+"	"+position[i].get(j).nucl+"	"+position[i].get(j).coverage);
					for (int k=0;k<position[i].get(j).samples.length;k++){
						output.write("	");
						if(position[i].get(j).samples[k].size()>0){
							output.write(""+position[i].get(j).samples[k].get(0));
							if(position[i].get(j).samples[k].size()>1){
								for (int l=1;l<position[i].get(j).samples[k].size();l++){
									output.write(";"+position[i].get(j).samples[k].get(l));
								}
							}
						}	
					}
					
					if(position[i].get(j).as_no!=-1){
						output.write("	"+position[i].get(j).as_ref+"	"+position[i].get(j).as_no+"	"+position[i].get(j).as_tumor1+"	"+position[i].get(j).as_tumor2+"	"+position[i].get(j).as_tumor3);
					}
					output.newLine();
					
				}
				output.close();
			}
			
			
		}
		catch(Exception e){
			StackTraceElement[] aa=e.getStackTrace();
			for (int i=0;i<aa.length;i++){
				System.out.println(i+"	"+aa[i].getLineNumber());
			}
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
	
	public static boolean isNucleotide(String s){
		if(s.equals("A")){
			return true;
		}
		else if(s.equals("C")){
			return true;
		}
		else if(s.equals("G")){
			return true;
		}
		else if(s.equals("T")){
			return true;
		}
		return false;
	}
	
	public static int type(String ref, String tumor){//conversion to the type index
		
		if(ref.equals("A")){
			if(tumor.equals("C")){
				return 1;
			}
			else if(tumor.equals("G")){
				return 2;
			}
			else if(tumor.equals("T")){
				return 0;
			}
		}
		else if(ref.equals("C")){
			if(tumor.equals("A")){
				return 0;
			}
			else if(tumor.equals("G")){
				return 1;
			}
			else if(tumor.equals("T")){
				return 2;
			}
		}
		else if(ref.equals("G")){
			if(tumor.equals("A")){
				return 2;
			}
			else if(tumor.equals("C")){
				return 1;
			}
			else if(tumor.equals("T")){
				return 0;
			}
		}
		else if(ref.equals("T")){
			if(tumor.equals("A")){
				return 0;
			}
			else if(tumor.equals("C")){
				return 2;
			}
			else if(tumor.equals("G")){
				return 1;
			}
		}
		//System.out.println(ref+">"+tumor);
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
	
	
	//This thread reads for 1 chromosome the reference sequence file
	private static class Subthread extends Thread{
		String chr="";
		volatile boolean done=false;
		Hashtable<Integer,Integer> table=new Hashtable<Integer,Integer>();
		ArrayList<Position> position=new ArrayList<Position>();
		
		
		public void run(){
			
			done=false;
			try{
				FileInputStream in=new FileInputStream(file_peptide+chr+".txt");
				DataInputStream inn=new DataInputStream(in);
				BufferedReader input= new BufferedReader(new InputStreamReader(inn));
				
				FileInputStream in2=new FileInputStream(file_coverage+chr+".txt");
				DataInputStream inn2=new DataInputStream(in2);
				BufferedReader input2= new BufferedReader(new InputStreamReader(inn2));
				
				String s1="";
				String s2="";
				
				int kk=0;
				while((s1=input.readLine())!=null){
					
					s2=input2.readLine();
					String[] t1=s1.split("	");
					String[] t2=s2.split("	");
					
					if(!t1[0].equals(t2[0])){
						System.exit(0);
					}
					int pos=Integer.parseInt(t1[0]);
					String nucl=t1[1];
					double coverage=Double.parseDouble(t2[2]);
					if(t1.length<3){
						table.put(pos,kk);
						position.add(new Position(pos,nucl,coverage));
					}
					else{
						table.put(pos,kk);
						String as_ref=t1[2];
						int as_no=Integer.parseInt(t1[3]);
						String as_tumor1=t1[4];
						String as_tumor2=t1[5];
						String as_tumor3=t1[6];
						position.add(new Position(pos,nucl,coverage,as_ref,as_no,as_tumor1,as_tumor2,as_tumor3));
					}
					kk++;
				}
				input.close();
				input2.close();
				
			}
			catch(Exception e){
				StackTraceElement[] aa=e.getStackTrace();
				for (int i=0;i<aa.length;i++){
					System.out.println(i+"	"+aa[i].getLineNumber());
				}
				System.out.println(e);
			}
			done=true;
		}
	}
	
	//the object encodes the annotation of each position and the indices of the samples which have a mutation
	private static class Position{
		int pos=-1;
		String nucl="";
		double coverage=0;
		String as_ref="";
		int as_no=-1;
		String as_tumor1="";
		String as_tumor2="";
		String as_tumor3="";
		ArrayList<Integer>[] samples=new ArrayList[3];
		
		public Position(int pos, String nucl, double coverage){
			this.pos=pos;
			this.nucl=nucl;
			this.coverage=coverage;
			samples[0]=new ArrayList<Integer>();
			samples[1]=new ArrayList<Integer>();
			samples[2]=new ArrayList<Integer>();
		}
		public Position(int pos, String nucl, double coverage,String as_ref, int as_no,String as_tumor1, String as_tumor2, String as_tumor3){
			this.pos=pos;
			this.nucl=nucl;
			this.coverage=coverage;
			this.as_ref=as_ref;
			this.as_no=as_no;
			this.as_tumor1=as_tumor1;
			this.as_tumor2=as_tumor2;
			this.as_tumor3=as_tumor3;
			samples[0]=new ArrayList<Integer>();
			samples[1]=new ArrayList<Integer>();
			samples[2]=new ArrayList<Integer>();
		}
		public void add(int sample_index,int type){
			samples[type].add(sample_index);
		}
	}
}
