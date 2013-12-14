(ns loom.test.clustering
  (:use [loom graph clustering] :reload)
  (:use [clojure.test]))

(deftest test-partition-conversions
  (let [set-vector [#{:a :b :c} #{:d :e} #{:f :g :h :i} #{:j} #{:k}]
        node-map {:a 0 :b 0 :c 0, :d 1 :e 1, :f 2 :g 2 :h 2 :i 2, :j 3, :k 4}]
    (testing "Partition conversions"
      (are [expected got] (= expected got)
        []         (partition-map->vector {})
        {}         (partition-vector->map [])
        [#{:a}]    (partition-map->vector {:a 0})
        {:a 0}     (partition-vector->map [#{:a}])
        set-vector (partition-map->vector node-map)
        node-map   (partition-vector->map set-vector)))))
      
(deftest test-modularity
  (let [g (graph [:a :b] [:a :c] [:b :c] [:b :d] [:c :e] [:c :f]
                 [:d :e] [:d :f] [:d :g] [:e :g] [:f :g] [:g :h])]
    (testing "Modularity"
      (are [expected got] (= expected got)
        0     (modularity g (partition-vector->map [#{:a :b :c :d :e :f :g :h}]))
        14/64 (modularity g (partition-vector->map [#{:a :b :c} #{:d :e :f :g :h}]))
        (- (/ (apply + (map #(* % %) [2 3 4 4 3 3 4 1]))
              (* 24 24)))  (modularity g {:a 0 :b 1 :c 2 :d 3 :e 4 :f 5 :g 6 :h 7})))
    (testing "Fast Newman"
      (is (= [#{:a :b :c} #{:d :e :f :g :h}] (fast-newman g))))))