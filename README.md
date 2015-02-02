---------------------------
Notes by XuSong
---------------------------


## 1. estimation ##

运行参数: lda est [alpha] [k] [settings] [data] [random/seeded/*] [directory]，具体的:

	KDDcup13:     est 0.2 80 settings.txt KDDcup13/KDD_bow random output_data
	20newsgroups: est 0.01 50 settings.txt 20_newsgroup/stopword2/train_for_lda_150.data random output_data
	douban:       est 0.01 50 settings.txt douban/dim_r.txt random output_data

模型输出:

* .[model] :alpha   beta 
* .theta(也就是gamma)，gamma虽然是dir参数，但是k维gamma的相对大小能够反映k维theta的相对大小,theta可以用gamma的归一化来表示)
* .beta, sum(exp(a))=1（已验证）
* likehood，也就是P(W|alpha,beta)
       

## 2. inference ##

已知alpha，beta求theta，z

运行参数: lda inf [settings] [model] [data] [name]

	比如 inf settings.txt final ap_1-500_for_infer.dat name_infe

* 载入的[model]有：final.other(包括alpha)  final.beta
* 输出：name_infe-gamma.dat  name_infe-lda-lhood.dat
* theta(也就是gamma)   likehood（也就是P(W|alpha,beta)）

* 实例参考 my_run2_alpha=0.01_k=20/inf，经过比较可以发现，final.gamma和inf出的gamma拟合程度很好



原始数据分析：
ap.dat中的第一篇文档，对应ap.txt中的<DOCNO> AP901128-0100 </DOCNO> 


------------------------------------------------------------------------

一.统计结果

这篇文档中的高频词汇有12:7，149:7，1870:7，3888:4，1971:4，也就是police，school，teacher分别出现了7次，classroom，boys分别出现4次。

二.LDA结果

经过LDA的gamma中选择出最大的两个topic是19,20,这两个topic下的police、school确实比重较大，证实了计算结果。（参考my_run2_alpha=0.01_k=20）

三.一点体会

这一方面也体现了tfidf与LDA的最大区别，tfidf只是找到了这个单词在文档中的权重，并未体现这个word的类别属性，也就是topic属性。比如只是tfidf能够说明baidu这个单词的比重比较大，但是并为体现baidu这个单词代表什么样的语义。