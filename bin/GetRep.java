
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.HashSet;

class FilterSequence {
	//Paired-end Fastq
	public void filterPairEndFastqWithRepeatKmer(HashSet<String> repeatSet,String filePath1,String filePath2,String outFile1,String outFile2 ,double threshold,String mode) throws IOException{
		System.out.println("The size of repeatSet is : " + repeatSet.size());
		File file1 = new File(filePath1);
		File file2 = new File(filePath2);
		if(file1.isFile() && file1.exists()
			&& file2.isFile() && file2.exists()){
			InputStreamReader reader1 = new InputStreamReader(new FileInputStream(file1),"utf-8");
			BufferedReader bufferedReader1 = new BufferedReader(reader1);
			InputStreamReader reader2 = new InputStreamReader(new FileInputStream(file2),"utf-8");
			BufferedReader bufferedReader2 = new BufferedReader(reader2);
			
			int lineNum = 0;
			String line1;
			String line2;
			String temp1 = new String();
			String temp2 = new String();
			int readStatus1 = 0;
			int readStatus2 = 0;
			FileWriter writer1 = new FileWriter(outFile1,true);
			FileWriter writer2 = new FileWriter(outFile2,true);
			while((line1 = bufferedReader1.readLine()) != null && (line2 = bufferedReader2.readLine()) != null){
				lineNum++;
				if(line1.length() == 0 || line2.length() == 0){
					continue;
				}
				if(lineNum % 4 == 1){
					temp1 = line1;
					temp2 = line2;
				}else if(lineNum % 4 == 2){
					readStatus1 = judgeRead3(line1, repeatSet, threshold);
					readStatus2 = judgeRead3(line2, repeatSet, threshold);
					if(mode.equals("both")){
						if(readStatus1 == 1 && readStatus2 == 1){
							writer1.write(temp1 + "\n" + line1 + "\n");
							writer2.write(temp2 + "\n" + line2 + "\n");
						}
					}else if(mode.equals("either")){
						if(readStatus1 == 1 || readStatus2 == 1){
							writer1.write(temp1 + "\n" + line1 + "\n");
							writer2.write(temp2 + "\n" + line2 + "\n");
						}
					}else if(mode.equals("N")){
						if(readStatus1 == 1 && readStatus2 == 1){
							writer1.write(temp1 + "\n" + line1 + "\n");
							writer2.write(temp2 + "\n" + line2 + "\n");
						}else if(readStatus1 == 1){
							writer1.write(temp1 + "\n" + line1 + "\n");
							line2 = line2.replaceAll("[a-zA-Z]", "N");
							writer2.write(temp2 + "\n" + line2 + "\n");
						}else if(readStatus2 == 1){
							line1 = line1.replaceAll("[a-zA-Z]", "N");
							writer1.write(temp1 + "\n" + line1 + "\n");
							writer2.write(temp2 + "\n" + line2 + "\n");
						}
					}
					
				}else{
					if(mode.equals("both")){
						if(readStatus1 == 1 && readStatus2 == 1){
							writer1.write(line1 + "\n");
							writer2.write(line2 + "\n");
						}
					}else if(mode.equals("either")){
						if(readStatus1 == 1 || readStatus2 == 1){
							writer1.write(line1 + "\n");
							writer2.write(line2 + "\n");
						}
					}else if(mode.equals("N")){
						if(lineNum % 4 == 3 && (readStatus1 == 1 || readStatus2 == 1)){
							writer1.write(line1 + "\n");
							writer2.write(line2 + "\n");
						}
						
						if(lineNum % 4 == 0){
							if(readStatus1 == 1 && readStatus2 == 1){
								writer1.write(line1 + "\n");
								writer2.write(line2 + "\n");
							}else if(readStatus1 == 1){
								writer1.write(line1 + "\n");
								line2 = line2.replaceAll(".", "!");
								writer2.write(line2 + "\n");
							}else if(readStatus2 == 1){
								line1 = line1.replaceAll(".", "!");
								writer1.write(line1 + "\n");
								writer2.write(line2 + "\n");
							}
						}
					}
					
				}
			}
			writer1.close();
			writer2.close();
			bufferedReader1.close();
			bufferedReader2.close();
		}
	}
	
	//PairedEnd Fasta
	public void filterPairEndFastaWithRepeatKmer(HashSet<String> repeatSet,String filePath1,String filePath2,String outFile1,String outFile2 ,double threshold,String mode) throws IOException{
		System.out.println("The size of repeatSet is : " + repeatSet.size());
		File file1 = new File(filePath1);
		File file2 = new File(filePath2);
		if(file1.isFile() && file1.exists()
			&& file2.isFile() && file2.exists()){
			InputStreamReader reader1 = new InputStreamReader(new FileInputStream(file1),"utf-8");
			BufferedReader bufferedReader1 = new BufferedReader(reader1);
			InputStreamReader reader2 = new InputStreamReader(new FileInputStream(file2),"utf-8");
			BufferedReader bufferedReader2 = new BufferedReader(reader2);
			
			String line1;    												//each line of the fastq file1
			String line2;													//each line of the fastq file2
			String comment1 = new String();									//keep current comment line	of file1
			String comment2 = new String();									//keep current comment line of file2			
			int readStatus1 = 0;											//keep the judge result of current file1 nucleotide sequences
			int readStatus2 = 0;											//keep the judge result of current file2 nucleotide sequences
			FileWriter writer1 = new FileWriter(outFile1,true);
			FileWriter writer2 = new FileWriter(outFile2,true);
			StringBuffer seq1 = new StringBuffer();
			StringBuffer seq2 = new StringBuffer();
			boolean firflag = true;
			while(true){
				while((line1 = bufferedReader1.readLine()) != null){
					if(line1.charAt(0) == '>'){
						comment1 = line1;
						break;
					}
					seq1.append(line1);
				}
				while((line2 = bufferedReader2.readLine()) != null){
					if(line2.charAt(0) == '>'){
						comment2 = line2;
						break;
					}
					seq1.append(line2);
				}
				if(seq1.length() != 0 && seq2.length() != 0){
					readStatus1 = judgeRead3(seq1.toString(), repeatSet, threshold);
					readStatus2 = judgeRead3(seq2.toString(), repeatSet, threshold);
					if(mode.equals("both")){
						if(readStatus1 == 1 && readStatus2 == 1){
							writer1.write(comment1 + "\n" + seq1.toString() + "\n");
							writer2.write(comment2 + "\n" + seq2.toString() + "\n");
						}
					}else if(mode.equals("N")){
						if(readStatus1 == 1 && readStatus2 == 1){
							writer1.write(comment1 + "\n" + seq1.toString() + "\n");
							writer2.write(comment2 + "\n" + seq2.toString() + "\n");
						}else if(readStatus1 == 1){
							writer1.write(comment1 + "\n" + seq1.toString() + "\n");
							writer2.write(comment2 + "\n" + seq2.toString().replaceAll("[a-zA-Z]", "N") + "\n");
						}else if(readStatus2 == 1){
							writer1.write(comment1 + "\n" + seq1.toString().replaceAll("[a-zA-Z]", "N") + "\n");
							writer2.write(comment2 + "\n" + seq2.toString() + "\n");
						}
					}
					seq1.setLength(0);
					seq2.setLength(0);
				}else if(firflag == true){
					break;
				}
				firflag = true;
			}
			
			writer1.close();
			writer2.close();
			bufferedReader1.close();
			bufferedReader2.close();
		}
	}
	
	//Single Fastq
	public void filterSingleFastqWithRepeatKmer(HashSet<String> repeatSet,String filePath,String outFile,double threshold) throws IOException{
		System.out.println("The size of repeatSet is : " + repeatSet.size());
		File file = new File(filePath);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
            BufferedReader bufferedReader = new BufferedReader(reader);
            String line;   														//each line of the fastq file
            int lineNum = 0;
            String comment = new String(); 										//keep current comment line
            int readStatus = 0;													//keep the judge result of current nucleotide sequences
            FileWriter writer = new FileWriter(outFile,true);
            while((line = bufferedReader.readLine()) != null){
            	lineNum++;
            	if(line.length() == 0){
            		continue;
            	}
            	
            	if(lineNum % 4 == 1){
            		comment = line;
            	}else if(lineNum % 4 == 2){
            		readStatus = judgeRead3(line, repeatSet, threshold);
            		if(readStatus == 1){
            			writer.write(comment + "\n" + line + "\n");
            		}
            	}else{
            		if(readStatus == 1){
            			writer.write(line + "\n");
            		}
            	}
            }
            writer.close();
            bufferedReader.close();
		}
	}
	
	//Single Fasta
	public void filterSingleFastaWithRepeatKmer(HashSet<String> repeatSet,String filePath,String outFile,double threshold) throws IOException{
		System.out.println("The size of repeatSet is : " + repeatSet.size());
		File file = new File(filePath);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
            BufferedReader bufferedReader = new BufferedReader(reader);
            String line;  														//each line of fasta file
            String comment = new String(); 										//keep current comment line
            int readStatus = 0; 												//keep the judge result of current nucleotide sequences
            FileWriter writer = new FileWriter(outFile,true);
            StringBuffer oneRead = new StringBuffer(); 							//keep current nucleotide sequences
            while((line = bufferedReader.readLine()) != null){
            	if(line.length() == 0){
            		continue;
            	}
            	
            	if(line.charAt(0) == '>' || line.charAt(0) == ';'){
            		if(oneRead.length() != 0){
            			readStatus = judgeRead3(oneRead.toString(), repeatSet, threshold);
            			if(readStatus == 1){
            				writer.write(comment + "\n" + oneRead.toString() + "\n");
            			}
            			oneRead.setLength(0);
            		}
            		comment = line;
            	}else{
            		oneRead.append(line);
            	}
            }
            
            //process the last nucleotide sequences
            if(oneRead.length() != 0 && comment.length() != 0){
    			readStatus = judgeRead3(oneRead.toString(), repeatSet, threshold);
    			if(readStatus == 1){
    				writer.write(comment + "\n" + oneRead.toString() + "\n");
    			}
    		}
            
            writer.close();
            bufferedReader.close();
		}
	}
	
	public void GeneratePairedEndLibraryWithRepeatKmer(
			HashSet<String> repeatSet,String pairedFile1,String pairedFile2,String outputPath,double threshold){
		
	}
	
	//all coverage
	public int judgeRead(String sequence, HashSet<String> repeatSet,double t){
		int kmerLength = 31;
		int length = sequence.length();
		double threshold = length * t;
		if(length < kmerLength){
			return 0;
		}
		int rlength = 0;
		int cindex = 0;
		for(int i = 0;i <= length - kmerLength;i++){
			if(repeatSet.contains(sequence.substring(i,i+kmerLength))){
				if(rlength == 0){
					rlength += kmerLength;
				}
				else if(i - cindex > kmerLength){
					rlength += kmerLength;
				}else{
					rlength += (i - cindex);
				}
				cindex = i;
				if(rlength >= threshold){
					return 1;
				}
			}
		}
		return 0;
	}
	
	//Nov 27
	public int judgeRead3(String sequence,HashSet<String> repeatSet,double t){
		int kmerLength = 31;
		int length = sequence.length();
		if(length < kmerLength){
			return 0;
		}
		double threshold = (length - kmerLength + 1) * t;
		int count = 0;
		if(!repeatSet.contains(sequence.substring(0, kmerLength)) || !repeatSet.contains(sequence.substring(length - kmerLength))){
			return 0;
		}
		for(int i = 0;i <= length - kmerLength;i++){
			if(repeatSet.contains(sequence.substring(i,i+kmerLength))){
				count++;
			}
			if(count >= threshold){
				//System.out.println(count);
				return 1;
			}
		}
		//System.out.println(count);
		return 0;
	}
	
	public int judgeRead2(String sequence,HashSet<String> repeatSet,double t){
		int kmerLength = 31;
		int length = sequence.length();
		if(length < kmerLength){
			return 0;
		}
		double threshold = (length - kmerLength + 1) * t;
		int count = 0;
		for(int i = 0;i <= length - kmerLength;i++){
			if(repeatSet.contains(sequence.substring(i,i+kmerLength))){
				count++;
			}
			if(count >= threshold){
				//System.out.println(count);
				return 1;
			}
		}
		//System.out.println(count);
		return 0;
	}
	
	public HashSet<String> getRepeatKmerSet(String repeatFilePath) throws IOException{
		File repeatFile = new File(repeatFilePath);
		HashSet<String> set = new HashSet<>();
		if(repeatFile.isFile() && repeatFile.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(repeatFile),"utf-8");
			BufferedReader bufferedReader = new BufferedReader(reader);
			String line;
			while((line = bufferedReader.readLine()) != null){
				if(line.length() == 0){
					continue;
				}
				if(line.charAt(0) != '>'){
					set.add(line);
					set.add(reverse(line));
				}
			}
			bufferedReader.close();
		}
		return set;
	}
	
	public String reverse(String s) {
		int length = s.length();
		StringBuilder reverseSb = new StringBuilder();
		for (int i = 0; i < length; i++){
			 if(s.charAt(i) == 'A')
			 {	
				 reverseSb.insert(0,"T");
			 }
			 else if(s.charAt(i) == 'T')
			 {
				 reverseSb.insert(0,"A");
			 }
			 else if(s.charAt(i) == 'G')
			 {
				 reverseSb.insert(0,"C");
			 }
			 else if(s.charAt(i) == 'C')
			 {
		         reverseSb.insert(0, "G");
			 }else{
				 reverseSb.insert(0, "N");
			 }
	    }
		return reverseSb.toString();
	}
	
	public static void main(String[] args) throws IOException {
		String singleOrPairEnd = args[0];// "S" or "P"
		String fastqOrfasta = args[1];// "fastq" or "fasta"
		FilterSequence fs = new FilterSequence();
		if(singleOrPairEnd.equals("S")){
			String fastq = args[2];
			String outFile = args[3];
			String repeatFilePath = args[4];
			double threshold = Double.parseDouble(args[5]);
			if(fastqOrfasta.equals("fastq")){
				fs.filterSingleFastqWithRepeatKmer(fs.getRepeatKmerSet(repeatFilePath), fastq, outFile, threshold);
			}else if(fastqOrfasta.equals("fasta")){
				fs.filterSingleFastaWithRepeatKmer(fs.getRepeatKmerSet(repeatFilePath), fastq, outFile, threshold);
			}
		}else if(singleOrPairEnd.equals("P")){
			String fastq1 = args[2];
			String fastq2 = args[3];
			String outFile1 = args[4];
			String outFile2 = args[5];
			String repeatFilePath = args[6];
			double threshold = Double.parseDouble(args[7]);
			String mode = args[8];// both,either,N
			if(fastqOrfasta.equals("fastq")){
				fs.filterPairEndFastqWithRepeatKmer(fs.getRepeatKmerSet(repeatFilePath),fastq1, fastq2, outFile1, outFile2,threshold,mode);
			}else if(fastqOrfasta.equals("fasta")){
				fs.filterPairEndFastaWithRepeatKmer(fs.getRepeatKmerSet(repeatFilePath),fastq1, fastq2, outFile1, outFile2,threshold,mode);
			}
		}
		
	}
}

class GetFreKmerWithThreshold{
	public static void getKmerFile(String kmerFile ,int th,String outputFile) throws IOException{
		File file = new File(kmerFile);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
			BufferedReader bReader = new BufferedReader(reader);
			String line = new String();
			boolean flag = false;
			FileWriter writer = new FileWriter(outputFile,true);
			while((line = bReader.readLine()) != null){
				if(line.charAt(0) == '>'){
					int count = Integer.parseInt(line.split(">")[1]);
					if(count >= th){
						flag = true;
					}else{
						flag = false;
					}
				}
				
				if(flag){
					writer.write(line + "\n");
				}
			}
			writer.close();
			bReader.close();
		}
	}
	
	public static void main(String[] args) throws IOException {
		String kmerFile = args[0];
		int th = Integer.parseInt(args[1]);
		String outputFile = args[2];
		getKmerFile(kmerFile,th,outputFile);
	}
}

class SplitContigsWithLengthTH{
	private static void SplitContigs(String contigFile,int threshold){
		File file = new File(contigFile);
		if(file.isFile() && file.exists()){
			InputStreamReader reader;
			try {
				reader = new InputStreamReader(new FileInputStream(file),"utf-8");
				BufferedReader bufferedReader = new BufferedReader(reader);
				String line = "";
				Boolean flag = false;
				FileWriter writer1 = new FileWriter("contigs_moreThan"+ threshold + ".fasta",true);
				FileWriter writer2 = new FileWriter("contigs_lessThan"+ threshold + ".fasta",true);
				while((line=bufferedReader.readLine()) != null){
					if(line.charAt(0) == '>'){
						if(Integer.parseInt(line.split("_")[3]) > threshold){
							flag = true;
						}else{
							flag = false;
						}
					}

					if(flag == true){
						writer1.write(line + "\n");
					}else{
						writer2.write(line + "\n");
					}
				}
				bufferedReader.close();
				writer1.close();
				writer2.close();
			} catch (Exception e) {
				System.out.println("SplitContigs Error");
			}
			
		}
	}
	
	/*public static String parseDirectory(String filePath){
		String[] strs = filePath.split("/");
		StringBuffer sb = new StringBuffer();
		for(int i = 0;i < strs.length - 1;i++){
			sb.append(str[i]);
		}
	}*/
	
	public static void main(String[] args) {
		String contigFile = args[0];
		int threshold = Integer.parseInt(args[1]);
		SplitContigs(contigFile, threshold);
	}
}

class FilterContigs{
	
	public static void main(String[] args) throws IOException {
		String contigsFilePath = args[0];
		String afterFilterFilePath = args[1];
		String lengthThreshold = args[2];
		String coverageThreshold = args[3];
		doFilterContigs(contigsFilePath, afterFilterFilePath, Integer.parseInt(lengthThreshold), Double.parseDouble(coverageThreshold));
	}
	
	public static void doFilterContigs(String contigsFilePath, String afterFilterFilePath, int lengthThreshold, double coverageThreshold) throws IOException{
		File file = new File(contigsFilePath);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
			BufferedReader bufferReader = new BufferedReader(reader);
			String currentLine;
			boolean isReserve = false;
			FileWriter writer = new FileWriter(afterFilterFilePath,true);
			while((currentLine = bufferReader.readLine()) != null){
				if(currentLine.charAt(0) == '>'){
					String[] strs = currentLine.split("_");
					if(Integer.parseInt(strs[3]) < lengthThreshold || Double.parseDouble(strs[5]) < coverageThreshold){
						isReserve = false;
					}else{
						writer.write(currentLine + "\n");
						isReserve = true;
					}
				}else{
					if(isReserve == true){
						writer.write(currentLine + "\n");
					}
				}
			}
			writer.close();
		}
	}
	
}


class FindThresholdFromHisto{
	public static void main(String[] args) throws IOException {
		String histoFile = args[0];
		float aveReadsLen = Float.parseFloat(args[1]);
		String thresholdFile = args[2];
		findMaxPoint(histoFile, aveReadsLen, thresholdFile);
	}
	
	public static int findMaxPoint(String histoFile, float aveReadsLen, String thresholdFile) throws IOException{
		File file = new File(histoFile);
		int[] x = new int[1000];
		int[] y = new int[1000];
		int[] diff = new int[1000];
		int maxPos = 0;
		int maxFre = 0;
		int lineNum = 0;
		int beginPoint = 0;
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
			BufferedReader bufferReader = new BufferedReader(reader);
			String currentLine;
			while((currentLine = bufferReader.readLine()) != null && lineNum < 1000){
				x[lineNum] = Integer.parseInt(currentLine.split(" ")[0]);
				y[lineNum] = Integer.parseInt(currentLine.split(" ")[1]);
				lineNum++;
			}
			bufferReader.close();
		}
		
		for(int i = 5;i < 1000;i++){
			if(y[i] - y[i-1] > 0){
				beginPoint = i - 1;
				break;
			}
		}
		
		if(beginPoint == 0){
			return -1;
		}
		
		for(int i = beginPoint;i < 1000;i++){
			if(y[i] > maxFre){
				maxFre = y[i];
				maxPos = x[i];
			}
			diff[i] = y[i] - y[i - 1];
		}
		
//		System.out.println("Max position is " + maxPos);
		
		int maxDiff = 0;
		int maxDiffPos = 0;
		for(int i = beginPoint + 3;i < 1000 - 4;i++){
			//System.out.print(i);
			int left = diff[i - 3] + diff[i - 2] + diff[i - 1] + diff[i];
			int right = diff[i + 1] + diff[i + 2] + diff[i + 3] +diff[i + 4];
			if(left > 0 && right < 0){
				if(left - right > maxDiff){
					maxDiff = left - right;
					maxDiffPos = i;
				}
			}
		}
		
//		System.out.println("Max diff position is " + maxDiffPos);
		if(Math.abs(maxPos - maxDiffPos) <= 10){
			int threshold = Math.round((float)(maxPos * aveReadsLen) / (aveReadsLen - 15 + 1));
			FileWriter writer = new FileWriter(thresholdFile,true);
			writer.write(String.valueOf(threshold) + "\n");
			writer.close();
			return maxPos;
		}
		return -1;
	}
	
}


class FilterContigsOfSoapdenovo{
	
	public static void main(String[] args) throws IOException {
		String contigsFilePath = args[0];
		String afterFilterFilePath = args[1];
		String lengthThreshold = args[2];
		//String coverageThreshold = args[3];
		doFilterContigs(contigsFilePath, afterFilterFilePath, Integer.parseInt(lengthThreshold));
	}
	
	public static void doFilterContigs(String contigsFilePath, String afterFilterFilePath, int lengthThreshold) throws IOException{
		File file = new File(contigsFilePath);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
			BufferedReader bufferReader = new BufferedReader(reader);
			String currentLine;
			boolean isReserve = false;
			FileWriter writer = new FileWriter(afterFilterFilePath,true);
			while((currentLine = bufferReader.readLine()) != null){
				if(currentLine.charAt(0) == '>'){
					String[] strs = currentLine.split("\\s+");
					if(Integer.parseInt(strs[2]) < lengthThreshold){
						isReserve = false;
					}else{
						writer.write(currentLine + "\n");
						isReserve = true;
					}
				}else{
					if(isReserve == true){
						writer.write(currentLine + "\n");
					}
				}
			}
			writer.close();
		}
	}
	
}

class SplitContigByLength{
	public static void main(String[] args) throws IOException {
		String contigsFilePath = args[0];
		String longFilePath = args[1];
		String shortFilePath = args[2];
		int lengthThreshold = Integer.parseInt(args[3]);
		doFilterContigs(contigsFilePath, longFilePath, shortFilePath, lengthThreshold);
	}
	
	public static void doFilterContigs(String contigsFilePath, String longFilePath, String shortFilePath, int lengthThreshold) throws IOException{
		File file = new File(contigsFilePath);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
			BufferedReader bufferReader = new BufferedReader(reader);
			String currentLine;
			boolean isLong = false;
			FileWriter writer1 = new FileWriter(longFilePath,true);
			FileWriter writer2 = new FileWriter(shortFilePath,true);
			while((currentLine = bufferReader.readLine()) != null){
				if(currentLine.charAt(0) == '>'){
					String[] strs = currentLine.split("_");
					if(Integer.parseInt(strs[3]) < lengthThreshold){
						isLong = false;
						writer2.write(currentLine + "\n");
					}else{
						isLong = true;
						writer1.write(currentLine + "\n");
					}
				}else{
					if(isLong){
						writer1.write(currentLine + "\n");
					}else{
						writer2.write(currentLine + "\n");
					}
				}
			}
			writer1.close();
			writer2.close();
			bufferReader.close();
		}
	}

}

class CutContigEnd{
	public static void main(String[] args) throws IOException {
		String contigFilePath = args[0];
		double percent = Double.parseDouble(args[1]);
		String outputFilePath = args[2];
		cutContigsBothEnd(contigFilePath,percent,outputFilePath);
	}

	public static void cutContigsBothEnd(String contigFilePath, double percent, String outputFilePath) throws IOException{
		File file = new File(contigFilePath);
		if(file.isFile() && file.exists()){
			InputStreamReader reader = new InputStreamReader(new FileInputStream(file),"utf-8");
			BufferedReader bufferReader = new BufferedReader(reader);
			String currentLine;
			FileWriter writer = new FileWriter(outputFilePath,true);
			ArrayList<String> eachLineList = new ArrayList<>();
			while((currentLine = bufferReader.readLine()) != null){
				eachLineList.add(currentLine);
			}
			bufferReader.close();
			int cutEndLength = 0;
			StringBuilder sb = new StringBuilder();
			for(int i = 0;i < eachLineList.size();i++){
				String line = eachLineList.get(i);
				if(line.charAt(0) == '>'){
					int length = Integer.parseInt(line.split("_")[3]);
					cutEndLength = (int)Math.round((double)length * (1 - percent));
					cutEndLength = cutEndLength % 2 == 0 ? cutEndLength / 2 : (cutEndLength - 1)/2;
					writer.write(line + "\n");
					sb.setLength(0);
				}else{
					sb.append(line);
					if(i + 1 < eachLineList.size() && eachLineList.get(i+1).charAt(0) == '>'){
						String str = sb.substring(cutEndLength);
						str = str.substring(0,str.length() - cutEndLength);
						writer.write(str + "\n");
						sb.setLength(0);
					}
				}
			}
			writer.close();
		}
	}
}
