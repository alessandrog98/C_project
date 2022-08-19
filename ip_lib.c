/*
 *ID GRUPPO : 38 
 *Alessandro Gasparello 878169
 *Alessandro Simonato 882339
*/

#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){	
	unsigned int i,j,z;
        
    ip_mat *creata = (ip_mat *)malloc(sizeof(ip_mat));
    if (!creata){
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1);
    }
    
    creata->stat = (stats *)malloc(sizeof(stats)*k);
    if (!(creata->stat)){
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1);
    }
	
    creata->h = h;
    creata->w = w;    
    creata->k = k;
    
	creata->data=(float***)malloc(sizeof(float**)*(h));
    if (!(creata->data)){
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1);
    }	

    for (i=0 ; i<h; i++){
   
		creata->data[i] = (float**)malloc(sizeof(float*)*(w));
        if (!(creata->data[i])){
            printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
            exit(1);
        }    	
    
        for (j=0 ; j<w ; j++){
	        
            creata->data[i][j] = (float*)malloc(sizeof(float)*k);
            if (!(creata->data[i][j])){
                printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
                exit(1);
            }   		
            for (z=0 ; z<k ; z++){				
        	    set_val(creata,i,j,z,v); 
	        }	        
        }	
	}
    
    

    return creata;
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat  *a, unsigned int i,unsigned int j,unsigned int k, float v){
	if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}

void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
	for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
		for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf ("%f",get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_free(ip_mat *a){
	unsigned int i,j;
    
    if (a){    
        for (i=0 ; i < a->h; i++){
	    
		    for (j=0 ; j < a->w; j++){
                	
                free(a->data[i][j]);
	        }	
	    }
        
        for (i=0 ; i < a->h; i++){
		    
            free(a->data[i]);    
	    }

        free(a->data);
	    free(a->stat);
	    free(a);
    }

}

void compute_stats(ip_mat * t){
    unsigned int i,j,v;
    float acc=0.0,min=0.0,max=0.0;

    if (!t) {printf("Errore parametro a NULL !!!\n"); exit(1);}    
    
    for (v=0; v < t->k ; v++){ 
        acc = 0.0;
        for (i=0; i < t->h ; i++){
            for (j=0 ; j < t->w ; j++){
                acc += get_val(t ,i ,j ,v );;
       	    }		
        }
        t->stat[v].mean = acc/(t->h*t->w);  
    }
	
	
    for (v=0 ;v < t->k  ;v++ ){ 
        min = get_val(t ,0 ,0 ,v );
        for (i=0 ;i < t->h ;i++ ){
            for (j=0 ;j < t->w ;j++ ){
                if (min > get_val(t ,i ,j ,v )){
                    min = get_val(t ,i ,j ,v );    
       	        }
            }		
        }
        t->stat[v].min = min;    
    }

    for (v=0 ;v < t->k ;v++ ){ 
        max = get_val(t ,0 ,0 ,v );
        for (i=0 ;i < t->h ;i++ ){
            for (j=0 ;j < t->w ;j++ ){
                if (max < get_val(t ,i ,j ,v )){
                   
                    max = get_val(t ,i ,j ,v );    
       	        }
            }		
        }
        t->stat[v].max = max;    
    }
}

ip_mat * ip_mat_copy(ip_mat * in){
	ip_mat *copy;
	unsigned int r,c,z;
    copy = NULL;    
    
    if (!in) {printf("Errore parametro a NULL !!!\n"); exit(1);}    
    
    copy = ip_mat_create( in->h, in->w, in->k, 0.0);
    
    if(copy) {    
        for (z=0 ; z<in->k ; z++){	
            for (r=0 ; r<in->h ; r++){
		        for (c=0 ; c<in->w ; c++){
			    
			     	copy->data[r][c][z] = in->data[r][c][z]; 
			    }	
		    }
	    }	
	    
    	compute_stats(copy);
    	return copy;
    }
    
    printf("Qualcosa è andato storto nella creazione della copia di ip_mat!!! :'(\n");
    exit(1);    
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    
    unsigned int dim_row;
	unsigned int dim_col;
	unsigned int r,c,z,i,j;
	ip_mat *sub;
    sub = NULL;	
    
    if (!t) {printf("Errore parametro a NULL !!!\n"); exit(1);}  
	    
    if (col_end <= t->w && row_end <= t->h){        
        dim_row = (row_end-row_start);
	    dim_col = (col_end-col_start);
	    
	    sub = ip_mat_create(dim_row,dim_col,t->k,0.0);
	    
        if (sub){	
            for (z=0 ; z<sub->k ; z++){	
                i=0;        
                for (r=row_start ; r<row_end ; r++){
                    j=0;		    
                    for (c=col_start ; c<col_end ; c++){
			        
			         	sub->data[i][j][z] = t->data[r][c][z]; 
                        j++;			
                    }	
                    i++;		
                }
	        }	
	    
	        return sub;	
        }
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1);
    }
    printf("Errore!!! Non posso creare una sotto_matrice se le dimensioni fuoriescono dalla matrice di partenza\n");
    exit(1);
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    unsigned int dim;
	unsigned int r,c,z;
	ip_mat *out;
	
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}
    if (!b) {printf("Errore parametro a NULL !!!\n"); exit(1);}
 
    if ((a->h == b->h && a->w == b->w) || (a->w == b->w && a->k == b->k) || (a->h == b->h && a->k == b->k)){   
    
        out = NULL;        
        
        switch(dimensione){
		    
		    case 0 :	dim = a->h + b->h; 
					    out = ip_mat_create ( dim, a->w , a->k , 0.0);
					    
                        if(out){
                            for (z=0 ; z<a->k ; z++){	
                				for (r=0 ; r<a->h ; r++){
		            				for (c=0 ; c<a->w ; c++){
		            					set_val(out ,r ,c ,z ,get_val(a ,r ,c ,z ));
							        }
						        }
					        }
					        for (z=0 ; z<b->k ; z++){	
                				for (r=a->h ; r<dim ;r++ ){
		            				for (c=0 ; c<b->w ; c++){
								        set_val(out ,r ,c ,z ,get_val(b ,r-(a->h) ,c ,z ));
							        }
						        }
					        }
                        }                        
                        else {
                            printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
                            exit(1);
                        }    					    
                        break;
		    
		    case 1 :	dim = a->w + b->w; 
					    out = ip_mat_create ( a->h, dim , a->k , 0.0);
					    
                        if(out){
                            for (z=0 ; z<a->k ; z++){	
                				for (r=0 ; r<a->h ; r++){
		            				for (c=0 ; c<a->w ; c++){
		            					set_val(out ,r ,c ,z ,get_val(a ,r ,c ,z ));
							        }
						        }
					        }
					        for (z=0 ; z<b->k ; z++){	
                				for (r=0; r<b->h ; r++){
		            				for (c=a->w; c<dim; c++){
								        set_val(out ,r ,c ,z ,get_val(b ,r ,c-(a->w) ,z ));
							        }
							        
						        }
					        }
                        }                        
                        else {
                            printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
                            exit(1);
                        }    	
					    break;
			    
		    case 2 :	dim = a->k + b->k; 
					    out = ip_mat_create ( a->h, a->w , dim , 0.0);
					    
                        if(out){                        
                            for (z=0 ; z<a->k ; z++){	
                				for (r=0 ; r<a->h ; r++){
		            				for (c=0 ; c<out->w ; c++){
		            					set_val(out ,r ,c ,z ,get_val(a ,r ,c ,z ));
							        }
						        }
					        }
					        for (z=a->k ; z<dim ; z++){	
                				for (r=0; r<b->h ; r++){
		            				for (c=0; c<b->w; c++){
								        set_val(out ,r ,c ,z ,get_val(b ,r ,c ,z-(a->k)));
							        }
						        }
					        }
					     }
                         else {
                            printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
                            exit(1);
                         }    	
                         break;
	    }
    
	return out;	
    }
    
    printf("Errore!!! le matrici devono avere due dimensioni uguali\n");
    exit(1);
}

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
	
    unsigned int r,c,z;
	ip_mat *sum;
	sum = NULL;
	    
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}
    if (!b) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    if ((a->h == b->h) && (a->w == b->w) && (a->k == b->k)){        
        sum = ip_mat_create(a->h, a->w,a->k, 0.0);
	    
        if (sum){
            for (z=0 ; z<sum->k ; z++){	
                for (r=0 ; r<sum->h ; r++){
		            for (c=0 ; c<sum->w ; c++){
			        
			         	set_val(sum ,r ,c ,z ,get_val(a ,r ,c ,z ) + get_val(b ,r ,c ,z ));
			        }	
		        }
	        }
	        
            compute_stats(sum);   
            return sum;
        }    
        
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1);    
    }
    
    printf("Errore!!! le matrici devono avere dimensioni uguali\n");
    exit(1);    
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
	
    
    unsigned int r,c,z;
	ip_mat *sub;
	sub = NULL;
	 
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}
    if (!b) {printf("Errore parametro a NULL !!!\n"); exit(1);}
 
    if ((a->h == b->h) && (a->w == b->w) && (a->k == b->k)){  
        sub = ip_mat_create(a->h, a->w,a->k, 1.0);
	    
        if (sub){
            for (z=0 ; z<sub->k ; z++){	
                for (r=0 ; r<sub->h ; r++){
		            for (c=0 ; c<sub->w ; c++){
			        
			         	set_val(sub ,r ,c ,z ,get_val(a ,r ,c ,z ) - get_val(b ,r ,c ,z ));
			        }	
		        }
	        }
            
            compute_stats(sub); 
            return sub;
        }
        
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1);  
    }

    printf("Errore!!! le matrici devono avere dimensioni uguali\n");
    exit(1);   
}

ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    ip_mat *mean;       
    unsigned int i,j,v;
    mean = NULL;
        
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}
    if (!b) {printf("Errore parametro a NULL !!!\n"); exit(1);}    

    if((a->h == b->h) && (a->w == b->w) && (a->k == b->k)){
        mean = ip_mat_create(a->h,a->w,a->k,0.0);
    
        if (mean) {
            for (v = 0; v < a->k; v++){
                for (i = 0; i < a->h; i++){            
                    for (j = 0; j < a->w; j++){
                        set_val (mean,i,j,v,((get_val(a,i,j,v)+get_val(b,i,j,v)) / 2));
                    }     
                }
            }
            
            compute_stats(mean); 
            return mean;
        }
        
        printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
        exit(1); 
    }
    
    printf("Errore!!! le matrici devono essere quadrate\n");
    exit(1); 
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
	unsigned int r,co,z;
	ip_mat *mul;
	mul = NULL;
	
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    mul = ip_mat_create(a->h, a->w,a->k, 0.0);

    if (mul) { 
        for (z=0 ; z<mul->k ; z++){	
            for (r=0 ; r<mul->h ; r++){
		        for (co=0 ; co<mul->w ; co++){
			    
			     	set_val(mul,r,co,z,get_val(a,r,co,z) * c);
			    }	
		    }
	    }	
        
        return mul;
    }
    
    printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
    exit(1);
}

ip_mat * ip_mat_add_scalar(ip_mat *a, float c){
	unsigned int r,co,z;
	ip_mat *add;
	add = NULL;
	
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    add = ip_mat_create(a->h, a->w, a->k, 0.0);
	
    if (add) {
        for (z=0 ; z<add->k ; z++){	
            for (r=0 ; r<add->h ; r++){
		        for (co=0 ; co<add->w ; co++){
			    
			     	set_val(add,r,co,z,get_val(a,r,co,z) + c);
			    }	
		    }
	    }

        return add;
    }
    
    printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
    exit(1);
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    ip_mat *out;
    unsigned int i,j,v;    
    unsigned int acc=0;    
    out = NULL;
    
    if (!in) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    out = ip_mat_create(in->h,in->w,in->k,0.0);
    
    if (out) {    
        for (i = 0; i < in->h; i++){
            for (j = 0; j < in->w; j++){
                acc=0;
                for ( v=0; v< in->k ; v++){
                    acc += get_val(in ,i ,j ,v );
                }
                acc = acc/(in->k);        
                for ( v=0; v< in->k ; v++){
                    set_val(out ,i ,j ,v ,acc );       
                            
                }
            }
        }

        return out;
    }
    
    printf("Qualcosa è andato storto nella creazione della ip_mat!!! :'(\n");
    exit(1);
}
ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }
    
    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    
    return b;
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    ip_mat *out;
    float mula,mulb;
    unsigned int v,i,j;    
    out = NULL;

    
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}
    if (!b) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    if ((a->h == b->h) && (a->w == b->w) && (a->k == b->k)){    

        mula = alpha;
        mulb = 1.0-alpha;

        out = ip_mat_create(a->h ,a->w ,a->k ,0.0 );
    
        if (out) {
            for ( v = 0; v< out->k ; v++){
                for (i = 0; i < out->h; i++){
                    for (j = 0; j < out->w; j++){
                        
                        set_val(out ,i ,j ,v ,(get_val(a ,i ,j ,v )*mula) + (get_val(b ,i ,j ,v )*mulb) );
                    }
                }
            }
        }            
        return out;
    }
    printf("Errore!!! le matrici devono avere le medesime dimensioni\n");
    exit(1);
}

void clamp(ip_mat * t, float low, float high){
    unsigned int i,j,v;
    
    if (!t) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    for (i = 0; i < t->h; i++){
        for (j = 0; j < t->w; j++){
            for ( v=0; v< t->k ; v++){
                if (get_val(t,i,j,v) > high){
                    set_val(t,i,j,v,high);
                }
                else if (get_val(t,i,j,v) < low){
                    set_val(t,i,j,v,low);
                }
            }
        }
    }
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    ip_mat *out;
    out = NULL;    

    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}   

    out = ip_mat_add_scalar(a,bright);
       
    return out;
}

void ip_mat_init_random(ip_mat * t, float mean, float var){
	unsigned int r,c,z;	

    if (!t) {printf("Errore parametro a NULL !!!\n"); exit(1);}    
 
    for (z=0 ; z < t->k ; z++){
        for (r=0 ; r < t->h ; r++){
		    for (c=0 ; c < t->w ; c++){
			    
			    set_val(t ,r ,c ,z , get_val(t ,r ,c ,z) + (get_normal_random(mean,var/2)));
		    }	
	    }
	}	
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    ip_mat *corr;    
    corr = NULL;        
    
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}    

    corr = ip_mat_copy(a);    

    ip_mat_init_random(corr,0.0,amount);    
    return corr;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    ip_mat *conv;
    ip_mat *sub_mat;
    ip_mat *padd;    
    unsigned int i,j,v,h,w;    
    unsigned int ir,ic;    
    float result = 0.0;     /*variabile di supporto per salvare i risultati del filtro applicato alle sotto_matrici*/
    
    if (!a) {printf("Errore parametro a NULL !!!\n"); exit(1);}
    if (!f) {printf("Errore parametro a NULL !!!\n"); exit(1);}

    padd = NULL;
    conv = NULL;
    
    h = f->h;
    w = f->w;
    
    conv = ip_mat_copy(a);
    
    padd = ip_mat_padding( a, h, w);
    
    for (v = 0; v < padd->k ; v++){
        h = f->h;
        for (i = 0; h <= padd->h; i++){
            w = f->w;         
            sub_mat = NULL;            
            for (j = 0; w <= padd->w; j++){
                sub_mat = ip_mat_subset(padd, i, h, j, w);     /*prendiamo sotto_matrici di dimensioni pari a quelle del filtro da applicare*/
                result=0.0;
                for (ir = 0; ir < f->h; ir++){
                    for (ic = 0; ic < f->w; ic++){
                        
                        if (f->k == 1){     /*controllo livelli del filtro*/
                             result += get_val(sub_mat,ir,ic,v) * get_val(f,ir,ic,1);                        
                        }                    
                        else {
                             result += get_val(sub_mat,ir,ic,v) * get_val(f,ir,ic,v);                   
                        }                    
                    }
                }             
                ip_mat_free(sub_mat);     /*liberiamo dalla memoria ogni sub_mat utilizzata nel passo precedente*/
                set_val(conv,i,j,v,result);     /*applichiamo ad ogni cella della ip_mat i valori dopo la convoluzione*/
                w++;
            }
            h++;
        }
    }
    ip_mat_free(padd);
    compute_stats(conv);   
    return conv;
}

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
    ip_mat *padd;
    unsigned int h,w;
    unsigned int i,j,v;

    if (!a) {printf("Errore parametro a NULL"); exit(1);}

    padd = NULL;
    
    h = (pad_h-1);
    w = (pad_w-1);
    padd = ip_mat_create(h+a->h,w+a->w,a->k,0.0);  /*creiamo una nuova ip_mat aggiungendo ad ogni lato un bordo in base al filtro da applicare*/
       
    if(padd) {    
        for ( v = 0; v< padd->k ; v++){
            for (i = 0; i < padd->h; i++){
                for (j = 0; j < padd->w; j++){
                    if (i<(h/2)+(a->h) && j<(w/2)+(a->w) && i>=(h/2) && j>=(w/2)){
                        set_val( padd, i, j, v,get_val( a, i-(h/2), j-(w/2), v));   /*ricopiamo i valori della vecchia matrice in quella nuova tenendo conto dei bordi aggiunti*/
                    }
                }
            }
        }
        return padd;
    }
    
    printf("Qualcosa è andato storto nella creazione della ip_mat con i bordi aggiuntivi!!! :'(\n");
    exit(1);
}

ip_mat * create_sharpen_filter(){
    ip_mat *filter;
    unsigned int i,j,v;
    filter = NULL;
    
    filter = ip_mat_create(3,3,3,0.0);
    
    if(filter){
        for (v = 0; v < 3; v++){
            for (i = 0; i < 3; i++){
                for (j = 0; j < 3; j++){
                    if ((j!=0 || i!=0) && (j!=0 || i!=2) && (j!=2 || i!=0) && (j!=2 || i!=2)){
                            
                        if (j==1 && i==1){
                            set_val (filter,j,i,v,5.0);
                        }
                        else 
                            set_val (filter,j,i,v,-1.0);
                    }
                }
            }
        }    

        return filter;
    }
    
    printf("Qualcosa è andato storto nella creazione del filtro richiesto!!! :'(\n");
    exit(1);
}

ip_mat * create_edge_filter(){
    ip_mat *filter;    
    unsigned int i,j,v;
    filter = NULL;

    filter = ip_mat_create(3,3,3,0.0);
    
    if(filter){
        for (v = 0; v < 3; v++){
            for (i = 0; i < 3; i++){
                for (j = 0; j < 3; j++){
                    
                    if (j==1 && i==1){
                        set_val (filter,j,i,v,8.0);
                    }
                    else {
                        set_val (filter,j,i,v,-1.0);
                    }        
                }
            }
        }
        
        return filter;
    }
    
    printf("Qualcosa è andato storto nella creazione del filtro richiesto!!! :'(\n");
    exit(1);
}

ip_mat * create_emboss_filter(){
    ip_mat *filter;    
    unsigned int i,j,v;
    filter = NULL;

    filter = ip_mat_create(3,3,3,0.0);
    
    if(filter){    
        for (v = 0; v < 3; v++){
            for (i = 0; i < 3; i++){
                for (j = 0; j < 3; j++){
                    if ((i!=0 || j!=2) && (i!=2 || j!=0)){
                            
                        if ((i==1 && j==1) || (i==1 && j==2) || (i==2 && j==1)){
                            set_val (filter,i,j,v,1.0);
                        }
                        if ((i==0 && j==1) || (i==1 && j==0)){
                            set_val (filter,i,j,v,-1.0);
                        }
                        if ((i==0 && j==0)){
                            set_val (filter,i,j,v,-2.0);
                        }        
                        if ((i==2 && j==2)){
                            set_val (filter,i,j,v,2.0);
                        }                
                    }
                }
            }
        }    
        
        return filter;
    }
    
    printf("Qualcosa è andato storto nella creazione del filtro richiesto!!! :'(\n");
    exit(1);
}

ip_mat * create_average_filter(unsigned int h, unsigned int w, unsigned int k){
    
    if (h >= 3 && w >= 3 && (k >= 1 || k <= 3) && h%2 != 0 && w%2 != 0){    
        ip_mat *filter;    
        float c;    
        filter = NULL;
        
        c=1.0/(h*w);
        
        filter = ip_mat_create(w,h,k,c);
        
        if(filter){
        
            return filter;
        }
        printf("Qualcosa è andato storto nella creazione del filtro richiesto!!! :'(\n");
        exit(1);    
    }
    printf("Le dimensioni del filtro inserite non sono consentite!!! \n");
    exit(1);
}

void rescale(ip_mat * t, float new_max){    

    
    float new_res;    
    unsigned int i,j,v;
    
    if (!t) {printf("Errore creazione ip_mat !!!\n"); exit(1);}
        
    compute_stats(t);    

    for (v = 0; v < t->k; v++){
        for (i = 0; i < t->h; i++){
            for (j = 0; j < t->w; j++){
                
                new_res = ((get_val (t,i,j,v)-(t->stat[v].min))/((t->stat[v].max)-(t->stat[v].min))) * new_max;         
                                
                set_val (t,i,j,v,new_res);
            }     
        }
    }     
}

ip_mat * create_gaussian_filter(unsigned int h, unsigned int w, unsigned int k, float sigma){
    
    if (h >= 3 && w >= 3 && k >= 1){
        ip_mat *filter;
        unsigned int i,j,v;
        float val,somma = 0.0,cx,cy,x,y;    
        filter = NULL;    
        filter = ip_mat_create(w,h,k,0.0);    
        
        if(filter){
            cx = (w-1)/2;
            cy = (h-1)/2;        
        
            for (v = 0; v < filter->k; v++){    
                somma = 0.0;
                for (i = 0; i < filter->h; i++){        
                                       
                    for (j = 0; j < filter->w; j++){
                        
                        y = i - cy;
                        x = j - cx;
                        val = 1.0/(2*PI*pow(sigma, 2))* exp(((pow(x, 2) + pow(y, 2)) / ((-1.0)*(2.0 * pow(sigma, 2)))));
                        somma += val;            
                        set_val(filter,i,j,v,val);
                    }  
                } 
                
                for (i = 0; i < filter->h; i++){        
                                       
                    for (j = 0; j < filter->w; j++){
                        
                        set_val(filter,i,j,v,(get_val(filter,i,j,v))/somma);
                    }
                }
                            
            }
            return filter;
        }
        printf("Qualcosa è andato storto nella creazione del filtro richiesto!!! :'(\n");
        exit(1);        
    }
    printf("Le dimensioni del filtro inserite non sono consentite!!! \n");
    exit(1);    
}

/*funzione extra per ruotare le immagini*/

ip_mat * ip_mat_rotate(ip_mat * in,unsigned int degree){ /* degree = 0 -> mirroring -- degree = 360 -> mat_copy */
    ip_mat *rotate;
	unsigned int r,c,z,j,i;

    if (!in) {printf("Errore parametro a NULL !!!\n"); exit(1);}    
        
    rotate = NULL;    
       
    switch (degree){        
        
        case 0:    rotate = ip_mat_create( in->h, in->w, in->k, 0.0);
                   if (!rotate) {printf("Errore creazione ip_mat !!!\n"); exit(1);}                    
                   for (z=0 ; z<in->k ; z++){	                        
                       for (r=0 ; r<in->h ; r++){
		                   j = (in->w) - 1;
                           for (c=0 ; c<in->w ; c++){
			            
			             	   set_val(rotate ,r ,c , z, get_val(in ,r ,j ,z )); 
                               j--;                                
                                			    
                            }		                                     
                        }
                    }
                    break;  
        
        case 90:    rotate = ip_mat_create( in->w, in->h, in->k, 0.0);
                    if (!rotate) {printf("Errore creazione ip_mat !!!\n"); exit(1);}                    
                    for (z=0 ; z<in->k ; z++){	
                        j = (in->h) - 1;                        
                        for (r=0 ; r<in->h ; r++){
		                    i = 0;
                            for (c=0 ; c<in->w ; c++){
			            
			             	    set_val(rotate ,i ,j , z, get_val(in ,r ,c ,z )); 
                                i++;                                
                                			    
                            }		                
                            j--;                        
                        }
                    }
                    break;	    
        
        case 180:   rotate = ip_mat_create( in->h, in->w, in->k, 0.0);
                    if (!rotate) {printf("Errore creazione ip_mat !!!\n"); exit(1);}                    
                    for (z=0 ; z<in->k ; z++){	
                        i = (in->h) - 1;                       
                        for (r=0 ; r<in->h ; r++){
                            for (c=0 ; c<in->w ; c++){
			            
			             	    set_val(rotate ,i ,c , z, get_val(in ,r ,c ,z )); 			    
                            }	
		                i--;
                        }
                    }
                    break;
	
	    
        case 270:   rotate = ip_mat_create( in->w, in->h, in->k, 0.0);
                    if (!rotate) {printf("Errore creazione ip_mat !!!\n"); exit(1);}                   
                    for (z=0 ; z<in->k ; z++){	
                        j = 0;                        
                        for (r=0 ; r<in->h ; r++){
		                    i = (in->w) - 1;
                            for (c=0 ; c<in->w ; c++){
			            
			             	    set_val(rotate ,i ,j , z, get_val(in ,r ,c ,z )); 
                                i--;                                
                                			    
                            }		                
                            j++;                        
                        }
                    }
                    break;
        
        case 360:   rotate = ip_mat_copy(in);
                    break;
    }	

    return rotate;
}    

