(ns ^{:doc "Graph-generating functions"
      :author "Justin Kramer"}
  loom.gen
  (:use [loom.graph :only [weighted? directed? add-nodes* add-edges* nodes edges]]
        [loom.metrics :only [degrees]]))

(defn gen-rand
  "Adds num-nodes nodes and approximately num-edges edges to graph g. Nodes
  used for each edge are chosen at random and may be chosen more than once."
  [g num-nodes num-edges & {:keys [min-weight max-weight loops seed]
                            :or {min-weight 1
                                 max-weight 1
                                 loops false
                                 seed (System/nanoTime)}}]
  (let [rnd (java.util.Random. seed)
        rand-w #(+ (.nextInt rnd (- max-weight min-weight)) min-weight)
        rand-n #(.nextInt rnd num-nodes)
        weighted? (weighted? g)
        nodes (range num-nodes)
        edges (for [_ (range num-edges)
                    :let [n1 (rand-n) n2 (rand-n)]
                    :when (or loops (not= n1 n2))]
                (if weighted?
                  [n1 n2 (rand-w)]
                  [n1 n2]))]
    (-> g
        (add-nodes* nodes)
        (add-edges* edges))))

(defn gen-rand-p
  "Adds num-nodes nodes to graph g with the probability p of an edge between
  each node."
  [g num-nodes p & {:keys [min-weight max-weight loops seed]
                    :or {min-weight 1
                         max-weight 1
                         loops false
                         seed (System/nanoTime)}}]
  (let [rnd (java.util.Random. seed)
        rand-w #(+ (.nextInt rnd (- max-weight min-weight)) min-weight)
        directed? (directed? g)
        weighted? (weighted? g)
        nodes (range num-nodes)
        edges (for [n1 nodes n2 nodes
                    :when (and (if directed?
                                 (or loops (not= n1 n2))
                                 (or (> n1 n2)
                                     (and loops (= n1 n2))))
                               (> p (.nextDouble rnd)))]
                (if weighted?
                  [n1 n2 (rand-w)]
                  [n1 n2]))]
    (-> g
        (add-nodes* nodes)
        (add-edges* edges))))

(defn gen-preferential
  "Adds num-nodes to graph g, each with k edges. New nodes connect preferentially to heavily connected nodes.
  This model is known as preferential attachment, scale-free and Barabási-Albert.

  This implementation uses a roulette wheel selection process, where each node receives as many tickets as their
  degree. For each inserted vertex, k neighbors are chosen without repetition to avoid parallel edges."
  [g num-nodes k & {:keys [seed]
                    :or {seed (System/nanoTime)}}]
  (when (or (directed? g) (weighted? g))
    (throw (IllegalArgumentException. "Graph is directed and/or weighted")))
  (let [rnd (java.util.Random. seed)
	n (count (nodes g))
	base (if (< n k)
		(add-edges* g (for [n1 (range n k), n2 (concat (nodes g) (range n k))
				    :when (not= n1 n2)]
				[n1 n2]))
		g)]
    (letfn [(get-nth-edge [hist x]
	      (let [[v deg] (first hist) ; if x=-1, return first v found
		    next-x (- x deg)]
		(if (pos? next-x)
		  (recur (rest hist) next-x)
		  v)))
	    (put-edges [g v]
	      (:g (reduce (fn [{:keys [g degs num-edges] :as state} _]
			    (let [rand-edge (dec (.nextInt rnd (inc num-edges))) ; if num-edges=0, rand-edge=-1
				  u (get-nth-edge degs rand-edge)]
			      (-> state
				(update-in [:degs] dissoc u)
				(update-in [:num-edges] - (get degs u))
				(assoc :g (add-edges* g [[u v]])))))
			  {:g g :num-edges (count (edges g)) :degs (degrees g)} ; TODO: update state at each new edge insertion
			  (range k))))]
      (let [[min-node max-node] (if (< n k)
				  [k (- num-nodes n)]
				  [n (+ num-nodes n)])]
	(reduce put-edges base (range min-node max-node))))))
