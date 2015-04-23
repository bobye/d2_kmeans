from numpy import loadtxt,transpose,histogram,argmax,argmin,zeros
from sklearn.cluster import KMeans
from scipy.spatial.distance import cdist

def get_kmeans_classification(ratio = 100, k=120):
    filename = 'usps_blankout' + str(ratio);
    train = loadtxt(filename + '_train.mat', delimiter=','); train=transpose(train);
    train_label = loadtxt(filename + '_train.label');
    test = loadtxt(filename + '_test.mat', delimiter=','); test=transpose(test);
    test_label = loadtxt(filename + '_test.label');
    km = KMeans(k).fit(train);

    train_size = len(train_label);
    test_size = len(test_label);

    center_label = zeros(k)
    for j in range(0, k):
        hist,bin_edges = histogram([train_label[i] for i in range(0, train_size) if km.labels_[i] == j], range(0,11)) 
        center_label[j] = argmax(hist)

    cd=cdist(km.cluster_centers_, test);
    pred_label = [center_label[i] for i in argmin(cd, 0)];
    error = sum([1 for i in range(0,test_size) if not pred_label[i] == test_label[i]]) / float(test_size)
    print ratio,k,error


def main():
    ratio_arr = [30,40,50,60,70,80,90,100]
    k_arr = [120, 180, 240, 300, 360]
    for ratio in ratio_arr:
        for k in k_arr:
            get_kmeans_classification(ratio, k)
        

if __name__ == "__main__":
    main()
