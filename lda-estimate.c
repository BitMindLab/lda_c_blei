// (C) Copyright 2004, David M. Blei (blei [at] cs [dot] cmu [dot] edu)

// This file is part of LDA-C.

// LDA-C is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your
// option) any later version.

// LDA-C is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
// for more details.

// You should have received a copy of the GNU General Public License
// along with this program; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA

#include "lda-estimate.h"

/*
 * perform inference on a document and update sufficient statistics
 *
 */

double doc_e_step(document* doc, double* gamma, double** phi,
                  lda_model* model, lda_suffstats* ss)
{
    double likelihood;
    int n, k;

    // posterior inference

    likelihood = lda_inference(doc, model, gamma, phi);

    // update sufficient statistics

    double gamma_sum = 0;
    for (k = 0; k < model->num_topics; k++)
    {
        gamma_sum += gamma[k];
        ss->alpha_suffstats += digamma(gamma[k]);
    }
    ss->alpha_suffstats -= model->num_topics * digamma(gamma_sum);

    for (n = 0; n < doc->length; n++)
    {
        for (k = 0; k < model->num_topics; k++)
        {
            ss->class_word[k][doc->words[n]] += doc->counts[n]*phi[n][k];
            ss->class_total[k] += doc->counts[n]*phi[n][k];
        }
    }

    ss->num_docs = ss->num_docs + 1;

    return(likelihood);
}



/*
 * writes the word assignments line for a document to a file
 *
 */

void write_word_assignment(FILE* f, document* doc, double** phi, lda_model* model)
{
    int n;

    fprintf(f, "%03d", doc->length);
    for (n = 0; n < doc->length; n++)
    {
        fprintf(f, " %04d:%02d",
                doc->words[n], argmax(phi[n], model->num_topics));
    }
    fprintf(f, "\n");
    fflush(f);
}


/*
 * saves the gamma parameters of the current dataset
 *
 */

void save_gamma(char* filename, double** gamma, int num_docs, int num_topics)
{
    FILE* fileptr;
    int d, k;
    fileptr = fopen(filename, "w");

    for (d = 0; d < num_docs; d++)
    {
	fprintf(fileptr, "%5.10f", gamma[d][0]);
	for (k = 1; k < num_topics; k++)
	{
	    fprintf(fileptr, " %5.10f", gamma[d][k]);
	}
	fprintf(fileptr, "\n");
    }
    fclose(fileptr);
}


/*
 * run_em
 *
 */

void run_em(char* start, char* directory, corpus* corpus)
{
	printf("running EM\n");

    int d, n;
    lda_model *model = NULL;
    double **var_gamma, **phi;

    // allocate variational parameters

    var_gamma = malloc(sizeof(double*)*(corpus->num_docs));
    for (d = 0; d < corpus->num_docs; d++)
	var_gamma[d] = malloc(sizeof(double) * NTOPICS);

    int max_length = max_corpus_length(corpus); //指的是所有文档中，word的最大数目
    phi = malloc(sizeof(double*)*max_length);
    for (n = 0; n < max_length; n++)
	phi[n] = malloc(sizeof(double) * NTOPICS);

    // initialize model

    char filename[100];

    lda_suffstats* ss = NULL;
    if (strcmp(start, "seeded")==0)
    {
        model = new_lda_model(corpus->num_terms, NTOPICS);
        //"seeded" initializes each topic to a distribution smoothed from a randomly                chosen document;


        ss = new_lda_suffstats(model);
        corpus_initialize_ss(ss, model, corpus);
        lda_mle(model, ss, 0);
        model->alpha = INITIAL_ALPHA;
    }
    // "Random" initializes each topic randomly;
    else if (strcmp(start, "random")==0)
    {
        model = new_lda_model(corpus->num_terms, NTOPICS);
        ss = new_lda_suffstats(model);
        random_initialize_ss(ss, model);
        lda_mle(model, ss, 0);
        model->alpha = INITIAL_ALPHA;

        printf("alpha is %f\n",model->alpha );

    }
    else
   // or, you can specify a model name to load a pre-existing model as the
   //initial model (this is useful to continue EM   from where it left off).
    {
        model = load_lda_model(start);
        ss = new_lda_suffstats(model);
    }


    sprintf(filename,"%s/000",directory);
    save_lda_model(model, filename);

    // run expectation maximization

    int i = 0;
    double likelihood, likelihood_old = 0, converged = 1;
    sprintf(filename, "%s/likelihood.dat", directory);
    FILE* likelihood_file = fopen(filename, "w");

    while (((converged < 0) || (converged > EM_CONVERGED) || (i <= 2)) && (i <= EM_MAX_ITER))
    {
        i++; printf("**** em iteration %d ****\n", i);
        likelihood = 0;
        zero_initialize_ss(ss, model);

        // e-step

        for (d = 0; d < corpus->num_docs; d++)
        {
            if ((d % 1000) == 0) printf("document %d\n",d);
            likelihood += doc_e_step(&(corpus->docs[d]),
                                     var_gamma[d],
                                     phi,
                                     model,
                                     ss);
        }

        // m-step

        lda_mle(model, ss, ESTIMATE_ALPHA);

        // check for convergence

        converged = (likelihood_old - likelihood) / (likelihood_old);
        if (converged < 0) VAR_MAX_ITER = VAR_MAX_ITER * 2;
        likelihood_old = likelihood;

        // output model and likelihood

        fprintf(likelihood_file, "%10.10f\t%5.5e\n", likelihood, converged);
        fflush(likelihood_file);
        if ((i % LAG) == 0)
        {
            sprintf(filename,"%s/%03d",directory, i);
            save_lda_model(model, filename);
            sprintf(filename,"%s/%03d.gamma",directory, i);
            save_gamma(filename, var_gamma, corpus->num_docs, model->num_topics);
        }
    }

    // output the final model

    sprintf(filename,"%s/final",directory);
    save_lda_model(model, filename);
    sprintf(filename,"%s/final.gamma",directory);
    save_gamma(filename, var_gamma, corpus->num_docs, model->num_topics);

    // output the word assignments (for visualization)

    sprintf(filename, "%s/word-assignments.dat", directory);
    FILE* w_asgn_file = fopen(filename, "w");
    for (d = 0; d < corpus->num_docs; d++)
    {
        if ((d % 100) == 0) printf("final e step document %d\n",d);
        likelihood += lda_inference(&(corpus->docs[d]), model, var_gamma[d], phi);
        write_word_assignment(w_asgn_file, &(corpus->docs[d]), phi, model);
    }
    fclose(w_asgn_file);
    fclose(likelihood_file);
}


/*
 * read settings.
 *
 */

void read_settings(char* filename)
{
    FILE* fileptr;
    char alpha_action[100];
    fileptr = fopen(filename, "r");
    fscanf(fileptr, "var max iter %d\n", &VAR_MAX_ITER);
    fscanf(fileptr, "var convergence %f\n", &VAR_CONVERGED);
    fscanf(fileptr, "em max iter %d\n", &EM_MAX_ITER);
    fscanf(fileptr, "em convergence %f\n", &EM_CONVERGED);
    fscanf(fileptr, "alpha %s", alpha_action);  //有两个值，estimate和fixed，详见readme
    if (strcmp(alpha_action, "fixed")==0)
    {
	ESTIMATE_ALPHA = 0;
    }
    else
    {
	ESTIMATE_ALPHA = 1;
    }
    fclose(fileptr);
}


/*
 * inference only
 *
 */

void infer(char* model_root, char* save, corpus* corpus)
{
    FILE* fileptr;   //载入模型model_root，infer后的参数存储在save里面
    char filename[100];
    int i, d, n;
    lda_model *model;
    double **var_gamma, likelihood, **phi;   //gamma的size是：N*K
    document* doc;

    model = load_lda_model(model_root);    //本代码的输入模型名称是final

    var_gamma = malloc(sizeof(double*)*(corpus->num_docs));  //初始化gammma
    for (i = 0; i < corpus->num_docs; i++)
	var_gamma[i] = malloc(sizeof(double)*model->num_topics);

    sprintf(filename, "%s-lda-lhood.dat", save);  //进行字符串连接，得到filename=“name_infe-lda-lhood.dat”
    fileptr = fopen(filename, "w");
    for (d = 0; d < corpus->num_docs; d++)
    {
		if (((d % 100) == 0) && (d>0)) printf("document %d\n",d);

		doc = &(corpus->docs[d]);
		phi = (double**) malloc(sizeof(double*) * doc->length);
		for (n = 0; n < doc->length; n++)
			phi[n] = (double*) malloc(sizeof(double) * model->num_topics);
		likelihood = lda_inference(doc, model, var_gamma[d], phi);   //---------------唯一核心------------

		fprintf(fileptr, "%5.5f\n", likelihood);  //写入变量likehood
    }
    fclose(fileptr);

    //写入变量gamma
    sprintf(filename, "%s-gamma.dat", save);
    save_gamma(filename, var_gamma, corpus->num_docs, model->num_topics);
}


/*
 * update sufficient statistics
 *
 */



/*
 * main
 *
 */
void mytest1()
{	double log10=2.3026;
	double log25=3.2189;
	double log35=log_sum(log10,log25);
	printf("你应该看到log35等于3.5553\n%f\n",log35);

}

/*tgamma----------- True gamma function. */
/* gamma------------Obsolete alias for `lgamma'.  */
void mytest2()
{	//n=4，则gamma=6，lgamma=1.791759
	double a=4;
	double d=tgamma(a);
	double b=gamma(a);
	//double c=log_gamma(a);
	printf("tgamma(4)=%f\nloggamma(4)=%f\n",d,b);

}

void mytest3()
{
	//make_directory("mytest3");
}

void mytest()
{
	mytest1();
	mytest2();
	mytest3();

}



int main(int argc, char* argv[])
{
    // (est / inf) alpha k settings data (random / seed/ model) (directory / out)

    corpus* corpus;

    long t1;
    (void) time(&t1);
    seedMT(t1);
    // seedMT(4357U);
    printf("1254");
    mytest();

    printf("4851");
    printf("        lda inf [settings] [model] [data] [name]\n");
    if (argc > 1)
    {
        if (strcmp(argv[1], "est")==0)
        {
            INITIAL_ALPHA = atof(argv[2]);
            NTOPICS = atoi(argv[3]);
            read_settings(argv[4]);
            //printf("2");
            corpus = read_data(argv[5]);
            //注意，数据中，第一个数字不是word总数，而是不同word的个数
            make_directory(argv[7]);
            run_em(argv[6], argv[7], corpus);
        }
        if (strcmp(argv[1], "inf")==0)
        {
        	//printf("3");
            read_settings(argv[2]);  //读取配置文件
            corpus = read_data(argv[4]);  //读取new_doc的数据
            infer(argv[3], argv[5], corpus);  //执行infer
        }
    }
    else
    {
        printf("usage : lda est [initial alpha] [k] [settings] [data] [random/seeded/*] [directory]\n");
        printf("        lda inf [settings] [model] [data] [name]\n");



    }
    return(0);
}
