(ns ^{:doc "Clustering methods in networks"
      :author "Bruno Kim Medeiros Cesar"}
  loom.clustering
  (:use [loom.graph :only [weighted? directed? nodes edges has-edge? successors out-degree]]))

(defn partition-map->vector
  "Converts a partition given as a {node part} map to a vector of sets.
  A partition identifier should be a number between in [0,p), where p
  is the number of partitions."
  [part]
  (let [p (inc (apply max -1 (vals part)))]
    (reduce (fn [set-vector [node part-index]]
              (update-in set-vector [part-index] conj node))
            (vec (repeat p #{}))
            part)))

(defn partition-vector->map
  "Converts a partition given as a vector of sets to a map of the form {node part}."
  [part]
  (reduce (fn [node-map [i vs]]
            (reduce #(assoc %1 %2 i) node-map vs))
          {}
          (map vector (range) part)))

(defn modularity
  "Calculates the modularity score of a given graph partition.
  A partition is given as a map of {node partition}."
  [g part]
  (let [p (inc (apply max -1 (vals part)))
        zeros (vec (repeat p (vec (repeat p 0))))
        ; It's necessary to compute m as well, because (edges g) doesnt guarantee
        ; if edges appear duplicated in a simple graph.
        [m mat] (reduce (fn [[m mat] [n1 n2]]
                          [(inc m) (update-in mat [(get part n1) (get part n2)] inc)])
                        [0 zeros]
                        (for [n1 (nodes g), n2 (successors g n1)] [n1 n2]))
        a (mapv #(/ (apply + %) m) mat)]
    (apply + (for [i (range p)]
               (let [eii (/ (get-in mat [i i]) m)
                     ai (get a i)]
                 (- eii (* ai ai)))))))

(defn- transpose [mat] (apply mapv vector mat))
(defn- remove-nth [v i] (vec (concat (subvec v 0 i) (subvec v (inc i)))))
(defn- remove-symm-row [mat i] (-> mat (remove-nth i) transpose (remove-nth i)))
(defn- sum-rows [mat i j]
  (vec (reduce (fn [m k]
                 (let [mik (get-in m [i k]), mjk (get-in m [j k])]
                   (assoc-in m [i k] (+ mik mjk)))) 
               mat 
               (range (count mat)))))

(defn fast-newman-seq
  "Returns lazy sequence of successive hierarchical clustering in the given graph by greedily
  maximizing modularity.
  
  Graph must be undirected and unweighted."
  [g]
  (when (or (directed? g) (weighted? g))
    (throw (IllegalArgumentException. "Graph is directed and/or weighted")))
  (let [vs (vec (nodes g))
        n (count vs)
        degs (mapv #(out-degree g %) vs)
        m (apply + degs)]
    (letfn [(most-connected-neighbor [ei i a]
              (let [ai (get a i)]
                (reduce #(max-key :dq %1 %2)
                        (map (fn [j eij aj]
                               {:pair [i j] 
                                :dq (if (= i j) 
                                      Double/NEGATIVE_INFINITY
                                      (* 2 (- eij (* ai aj))))})
                             (range) ei a))))
            (most-connected-pair [mat a]
              (let [vs (pmap most-connected-neighbor mat (range) (repeat a))]
                (reduce #(max-key :dq %1 %2) vs)))
            (join-clusters [{:keys [mat q tree]}]
              (let [a (mapv #(apply + %) mat)
                    {:keys [pair dq]} (most-connected-pair mat a)
                    [i j] (sort pair)]
                {:mat (-> mat 
                        (sum-rows i j) 
                        transpose 
                        (sum-rows i j) 
                        (remove-symm-row j) 
                        (assoc-in [i i] 
                          (+ (get-in mat [i i]) (* 2 (get-in mat [i j])) (get-in mat [j j]))))
                 :q (+ q dq)
                 :tree (-> tree (remove-nth j) (assoc i [(tree i) (tree j)]))}))]
      (take n 
        (iterate join-clusters
          {:mat  (vec (for [n1 vs] 
                   (vec (for [n2 vs] 
                     (if (has-edge? g n1 n2) 
                       (/ 1 m) 
                       0)))))
           :q    (- (/ (apply + (map #(* % %) degs))
                       (* m m)))
           :tree vs})))))

(defn fast-newman
  "Returns a graph partition that maximizes modularity according to the Fast Newman greedy algorithm.
  The partition is a vector of sets."
  [g]
  (vec (->> (fast-newman-seq g)
         (apply max-key :q )
         :tree
         (map (comp set flatten)))))
