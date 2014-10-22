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

void trim_read( Read * read, int qual);
void free_read(Read * read);
int main (int argc, char *argv[])
{
 	if (argc==1){
		printf("Usage: ./fastq_qual_trimmer <INPUT> -q 20 -l 1\n");
		return 0;
	}
	int opt = 0;
	int min_read_len = 1;
	int quality_threshold=20;
	char * fastq_name =strdup(argv[1]);
	while((opt = getopt(argc, argv, "q:l:")) != -1){
		switch(opt){
			case 'q':
			quality_threshold=atoi(optarg);
			//printf("Q:%d\n",quality_threshold);
			break;
			case 'l':
			min_read_len = atoi(optarg);
			break;
		}
	}
	//parse fastq file
	FILE * fp = fopen(fastq_name,"read");
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
			trim_read(&read,quality_threshold);
			int read_len=strlen(read.seq);
			if (read_len >= min_read_len){
				printf("%s\n%s\n%s\n%s\n",read.id,read.seq,read.plus,read.qual);
			}
			is_read_parsed=0;
		}
	}
	//free_read(&read);
	free(fastq_name);
	free(tmp);
	fclose(fp);
	return 0;

}
void trim_read( Read * read, int qual_thresh){
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
}
void free_read(Read * read){
	free(read->seq);
	free(read->qual);
	free(read->plus);
	free(read->id);
}
