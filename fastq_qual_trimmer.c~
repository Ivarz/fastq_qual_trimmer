#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
typedef struct Reads
{
	char * id;
	char * seq;
	char * plus;
	char * qual;
} Read;
void mean_qual_filter(Read * read, int qual_thresh);
void trim_read( Read * read, int qual);
void trim_window_read( Read * read, int qual_thresh, int window_size);
void free_read(Read * read);
int main (int argc, char *argv[])
{
	int print_help=0;
	int opt = 0;
	int min_read_len = 1;
	int quality_threshold = 20;
	int windowed = 0;
	int both = 0;
	int window_size = 5;
	int mean_q_thresh = 0;
	char * opt_fastq_name = NULL;
	while((opt = getopt(argc, argv, "i:hq:l:wbs:m:")) != -1){
		switch(opt){
			case 'i':
			opt_fastq_name=optarg;
			break;
			case 'q':
			quality_threshold=atoi(optarg);
			//printf("Q:%d\n",quality_threshold);
			break;
			case 'l':
			min_read_len = atoi(optarg);
			break;
			case 'w':
			windowed = 1;
			break;
			case 's':
			window_size = atoi(optarg);
			case 'b':
			both = 1;
			break;
			case 'm':
			mean_q_thresh = atoi(optarg);
			break;
			case 'h':
			print_help=1;
			break;
		}
	}
	//parse fastq file

	FILE * fp;
	if (opt_fastq_name!=NULL){
		fp = fopen(opt_fastq_name,"read");
	} else {
		fp = stdin;
	}
	if (print_help==1 || argc == 1){
		printf("Read quality trimming. Currently supports only fastq files. Output goes to STDOUT\n"
		"Usage:"
		"fastq_qual_trimmer -i <fastq_file> <OPTIONS>\n"
		"where options are:\n"
		"\t\t-h: print this message\n"
		"\t\t-i STRING: input file\n"
		"\t\t-q INT: trimming quality threshold [20]\n"
		"\t\t-m INT: read mean quality threshold [0]\n"
		"\t\t-l INT: minimal read length threshold [1]\n"
		"\t\t-w: use sliding window trimming\n"
		"\t\t-s INT: sliding window size [5]\n"
		"\t\t-b: use sliding window and after that - default trimming\n");
		return 0;
	}
	int character;
	int is_id=1;
	int is_seq=0;
	int is_plus=0;
	int is_qual=0;
	int is_read_parsed=0;
	char * tmp;
	Read read;
	read.id = malloc(sizeof(char));
	read.seq = malloc(sizeof(char));
	read.qual = malloc(sizeof(char));
	read.plus = malloc(3*sizeof(char));
	int idx = 0;
	while ((character=fgetc(fp))!=EOF){
		if (is_id==1 && character!='\n' && character!=EOF){
			//Parse fastq id string
			read.id[idx] = character;
			tmp = realloc(read.id, (idx+2)*sizeof(char));
			if (tmp !=NULL){
        		read.id=tmp;
      		} else {
        		free(read.id);
        		printf("Error allocating memory!\n");
      		}
			read.id[idx+1]='\0';
			++idx;
		}
		if (character == '\n' && is_id==1){
			is_id=0;
			is_seq=1;
			is_plus=0;
			is_qual=0;
			idx=0;
			continue;
		}
		//parse fastq seq
		if (is_seq==1 && character!='\n' && character!=EOF){
			read.seq[idx] = character;
			tmp = realloc(read.seq, (idx+2)*sizeof(char));
			if (tmp !=NULL){
			    read.seq=tmp;
			} else {
				free(read.seq);
				printf("Error allocating memory!\n");
			}
			read.seq[idx+1]='\0';
			++idx;
		}
		if (character == '\n' && is_seq==1){
			is_id=0;
			is_seq=0;
			is_plus=1;
			is_qual=0;
			idx=0;
			continue;
		}
		//Parse fastq +
		if (is_plus==1 && character!='\n' && character!=EOF){
			strcpy(read.plus, "+");
		}
		if (character == '\n' && is_plus==1){
			is_id=0;
			is_seq=0;
			is_plus=0;
			is_qual=1;
			idx=0;
			continue;
		}
		if (is_qual==1 && character!='\n' && character!=EOF){
            read.qual[idx] = character;
            tmp = realloc(read.qual, (idx+2)*sizeof(char));
            if (tmp !=NULL){
                read.qual=tmp;
            } else {
                free(read.qual);
                printf("Error allocating memory!\n");
            }
			read.qual[idx+1]='\0';
            ++idx;
        }
		if (character == '\n' && is_qual==1){
			is_id=1;
			is_seq=0;
			is_plus=0;
			is_qual=0;
			idx=0;
			is_read_parsed=1;
		}
		if (is_read_parsed == 1){
			if (both == 1){
				trim_window_read(&read,quality_threshold,window_size);
				trim_read(&read,quality_threshold);
			} else if (windowed == 1){
				trim_window_read(&read,quality_threshold,window_size);
			} else {
				trim_read(&read,quality_threshold);
			}
			mean_qual_filter(&read,mean_q_thresh);
			int read_len=strlen(read.seq);
			if (read_len >= min_read_len){
				printf("%s\n%s\n%s\n%s\n",read.id,read.seq,read.plus,read.qual);
			}
			is_read_parsed=0;
		}
	}
	//free_read(&read);
	free(tmp);
	fclose(fp);
	return 0;

}
void trim_read( Read * read, int qual_thresh)
{
	int start_idx=0;
	char * sequence;
	char * quality;
	int last_idx = strlen(read->qual);
	last_idx-=1;
	int end_idx=last_idx;
	int i;
	int stop_start_check=0;
	int stop_end_check=0;
	for (i = 0; i <=last_idx; ++i){
		if (read->qual[i] >= (qual_thresh+33) && stop_start_check==0){
			start_idx = i;
			stop_start_check=1;
		}
		if (read->qual[last_idx-i] >= (qual_thresh+33) && stop_end_check==0){
			end_idx = last_idx - i;
			stop_end_check=1;
		}
	}
	//If read is below qual threshold
	if(stop_end_check==0 && stop_start_check==0){
		read->seq[0]='\0';
		read->qual[0]='\0';
		return;
	}
	sequence = malloc((end_idx+2)*sizeof(char));
	quality = malloc((end_idx+2)*sizeof(char));
	strncpy(sequence, read->seq+start_idx, end_idx-start_idx+1);
	strncpy(quality, read->qual+start_idx, end_idx-start_idx+1);
	sequence[end_idx-start_idx+1]='\0';
	quality[end_idx-start_idx+1]='\0';
	
//	printf("%d\t%d\n%s\n%s\n%s\n%s\n",start_idx,end_idx,read->seq,read->qual,sequence,quality);
	strcpy(read->seq,sequence);
	strcpy(read->qual,quality);
	free(sequence);
	free(quality);
	return;
}
void trim_window_read( Read * read, int qual_thresh, int window_size){
	int start_idx=0;
	char * sequence;
	char * quality;
	int last_idx = strlen(read->qual);
	last_idx-=1;
	int end_idx=last_idx;
	int i;
	int j;
	int sum_of_q;
	int stop_start_check=0;
	int stop_end_check=0;
	//printf("%s\n",read->seq);
	if ((last_idx+1) < window_size){
		window_size=last_idx+1;
	}
	//printf("%d\n",window_size);
	for (i = 0; i <=(last_idx-window_size-1); ++i){
		sum_of_q=0;
		for (j=i; j<=i+window_size-1; ++j){
			sum_of_q=sum_of_q+read->qual[j];
		}
		//printf("%f\n",sum_of_q/(float)window_size);
		if (sum_of_q/(float)window_size >= (qual_thresh+33) && stop_start_check==0){
            start_idx = i+window_size-1;
            stop_start_check=1;
        }
        if (sum_of_q/(float)window_size < (qual_thresh+33) && stop_end_check==0 && stop_start_check==1){
            end_idx = i;
            stop_end_check=1;
        }
	}
   //If read is below qual threshold
    if(stop_end_check==0 && stop_start_check==0){
        read->seq[0]='\0';
        read->qual[0]='\0';
        return;
    }
    sequence = malloc((end_idx+2)*sizeof(char));
    quality = malloc((end_idx+2)*sizeof(char));
    strncpy(sequence, read->seq+start_idx, end_idx-start_idx+1);
    strncpy(quality, read->qual+start_idx, end_idx-start_idx+1);
    sequence[end_idx-start_idx+1]='\0';
    quality[end_idx-start_idx+1]='\0';
    
//  printf("%d\t%d\n%s\n%s\n%s\n%s\n",start_idx,end_idx,read->seq,read->qual,sequence,quality);
    strcpy(read->seq,sequence);
    strcpy(read->qual,quality);
    free(sequence);
    free(quality);
	
}
void mean_qual_filter(Read * read, int qual_thresh){
    int start_idx=0;
    char * sequence;
    char * quality;
    int seq_len = strlen(read->qual);
    int last_idx=seq_len - 1;
	int i;
	int total_qual = 0;
	for (i=0; i<=last_idx; ++i){
		total_qual+=read->qual[i];
	}
	if( total_qual/(float)seq_len < qual_thresh+33){
		sequence = malloc((2)*sizeof(char));
		quality = malloc((2)*sizeof(char));
		sequence[0]='\0';
		quality[0]='\0';
    	strcpy(read->seq,sequence);
    	strcpy(read->qual,quality);
    	free(sequence);
    	free(quality);	
	}
}

void free_read(Read * read){
	free(read->seq);
	free(read->qual);
	free(read->plus);
	free(read->id);
}
