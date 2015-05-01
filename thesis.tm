<TeXmacs|1.99.2>

<style|<tuple|article|algorithm>>

<\body>
  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  \;

  <doc-data|<doc-title|Structure from Motion System in Modern Browsers>>

  Pu Xu

  <abstract-data|<\abstract>
    3D scene reconstruction is a crucial aspect of computer vision. In this
    paper, we will design a structure from motion system which take photos
    took by regular pinehole cameras as input and runs inside a modern
    bowser.\ 

    Photos took by regular digital camera i.e. pinhold camera model are the
    most common and aboundant material for scene reconstruction. On the other
    hand, the most pervasive computation platform is the web browser, which
    is compatible with virtually every platform regardless the difference of
    operating system or hardware. 3D scene reconstruction system on the most
    pervasive platform - web browser, using the most common source material -
    photos, minimized the requirements for this important task, has great
    potential.

    Keywords: structure from m<name|>otion, 3D scene analysis, javascript,
  </abstract>>

  <new-page>

  <section|Introduction>

  \;

  The title ``Structure from Motion System in Modern Browser'' implies two
  key aspects, ``Structure from Motion'' and ``Modern Browser'', and on both
  sides it has great potential.\ 

  On one hand, structure from motion has always been a classic problem in
  machine vision. It's the problem of converting images or videos into 3D
  scene model, which is the exact reverse of CG(Computer Graphics). In resent
  years, structure from motion has became increasingly relavent because of
  the rising of technologies like VR (Virtual Reality), 3D printing and GIS
  (Geography Information System). Reconstruct 3D scene directly from regular
  photos i.e. follow pinehole camera model instead of specialized equipments
  can dramatically decreased the requirements for this important task, and
  open up the door to potentially unlimited resource for 3D reconstruction.

  On the other hand, modern browser as a platform is gaining its popularity.
  It has became increasingly powerful over the years both in functionality
  and performance. Several novel features of the modern browser has made
  implemeting an SFM system possible. First is the WebWorkers. Javascript in
  browser is mostly a single threaded language. With WebWorker, parallel
  computing with clent-side javascript became possible. Another one is WebGL,
  which is the openGL interface for the web. It made the rendering of 3D
  model practical (SVG and canvas can also render 3D model, but not as
  convenient). Last but not the least, Indexed DB is also a very important
  element. It implements a light weight NoSQL database, enables object
  storage and file storage inside the browser. For projects dealing with a
  large dataset and need efficient data persistent such as WebSFM, they are
  no longer confined in the limited memory space, makes IndexedDB an
  extremely crucial feature.

  In sum, the increasing need for 3D reconstruction in the industry and the
  maturity of web platform have made the SFM system in modern browser
  relavant and possible.

  <section|Overview>

  <big-figure|<image|framework.jpeg|700px|500px||>|Framework of WebSFM>

  The main process of SFM system happens in the second thread dedicated for
  the serial logic of the SFM pipeline. The pipline is based on multiple
  exsiting implementations, such as <cite|handheld><cite|modeltheworld>. And
  the sparse bundle adjustment part is influenced by SBA<cite|lma>. WebSFM
  has 3 phases:

  First is feature extraction. Image features are extracted using
  Lowe's<cite|sift> SIFT. It uses Difference of Guassian of the grayscale
  image to approximate Laplassian of Guassian of the image. Feature points
  are first found in the Laplacian of Guassian then filtered, finally binded
  with a descriptor. Features are scale/orientation invariant, which mean
  they can be found and matched in arbitray scale and orientation.

  Second is camera registration. We use the feature coorespondence between
  images to calibrate the cameras i.e. estimate the camera parameters. First
  we choose tracks i.e. globally unique and consistent feature matches. Then
  set up the initial pair of cameras. Afterwards, cameras are calibrated
  incrementally one after another, untill no reliable camera is left.

  Finally is stereopsis. After cameras are calibrated, we can use stereo
  vision to obtain dense point cloud. For each pair of cameras, rectify their
  relative pose into standard stereo configuration, then search feature match
  along the scan line which coincides epipolar line because of rectification.
  Then merge the two view stereo together to form MVS (Multi-View Stereo).

  Optionally, we can reconstruct surface or create mesh from the dense point
  cloud.

  <section|Mathematical Context>

  <subsection|View Geometry>

  <subsubsection|Coordinate System>

  Points in the SfM pipline have multiple cooridates according to different
  frames of reference.

  <\itemize-dot>
    <item><code*|Image Coordinate> -- <math|x>

    <item><code*|Camera Coordinate> -- <math|X<rsub|cam>>

    <item><code*|World Coordinate> -- <name|<math|X>>
  </itemize-dot>

  Transition between these coordiate systems are accomplished by the use of
  \ homogenous coordinate and transition matrices.

  <\itemize>
    <item><code*|Rotation Matrix(R) and Translation Vector(t)> --
    <math|X<rsub|cam>=<around*|[|R\<nocomma\>\<nocomma\>,t|]>X>

    <item><code*|Calibration Matrix(K)> -- <math|x =K X<rsub|cam>>\ 

    <item><code*|Camera/Projection Matrix(P)> -- <math|x=P X = K<around*|[|R
    ,t|]>X>
  </itemize>

  <subsubsection|View Constrains>

  Two-view constrains:

  <\itemize>
    <item><code*|Fundamental Matrix(F)> --
    <math|<with|math-display|true|x<rsup|T><rsub|1>F
    x<rsub|2>=0<with|math-display|false|>>>

    <item><code*|Essential Matrix(E)> -- \ <math|<with|math-display|true|X<rsup|T><rsub|cam1>E
    X<rsub|cam2>=0<with|math-display|false|>>>
  </itemize>

  Fundamental matrix is used for un-calibrated cameras, essential matrix is
  used for calibrated cameras, their relation is <math|F=K<rsup|T><rsub|1>E
  K<rsub|2>> .

  <subsubsection|Normalized Eight Point Algorithm (DLT)>

  Recover fundamental matrix of a two-view pair.\ 

  To be more specific, the algorithm used here is normalized
  8-point-algorithm. it first normalize the image coordinates into
  <math|<around*|[|0,1|]>\<times\><around*|[|0,1|]>>.

  <\equation*>
    x=N x=<matrix|<tformat|<table|<row|<cell|1/width>|<cell|0>|<cell|0>>|<row|<cell|0>|<cell|1/height>|<cell|0>>|<row|<cell|0>|<cell|0>|<cell|1>>>>><matrix|<tformat|<table|<row|<cell|x>>|<row|<cell|y>>|<row|<cell|1>>>>>
  </equation*>

  the epipolar constrain of fundamental matrix and a pair of matching points
  can be represented as:\ 

  <\equation*>
    x<rsup|T> F x<rprime|'>=<matrix|<tformat|<table|<row|<cell|x>|<cell|y>|<cell|1>>>>><matrix|<tformat|<table|<row|<cell|f<rsub|11>>|<cell|f<rsub|12>>|<cell|f<rsub|13>>>|<row|<cell|f<rsub|21>>|<cell|f<rsub|22>>|<cell|f<rsub|23>>>|<row|<cell|f<rsub|31>>|<cell|f<rsub|32>>|<cell|f<rsub|33>>>>>><matrix|<tformat|<table|<row|<cell|x<rprime|'>>>|<row|<cell|y<rprime|'>>>|<row|<cell|1>>>>>=0
  </equation*>

  it can also be written as:

  <\equation*>
    x<rprime|'> x f<rsub|11> + x<rprime|'> y f<rsub|12> + x<rprime|'>
    f<rsub|13> + y<rprime|'> x f<rsub|21> + y<rprime|'> y f<rsub|22> +
    y<rprime|'> f<rsub|23> + x f<rsub|31> + y f<rsub|32> + f<rsub|33> = 0
  </equation*>

  or in matrix form:

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|x<rprime|'> x>|<cell|x<rprime|'>
    y>|<cell| x<rprime|'>>|<cell|y<rprime|'> x >|<cell| y<rprime|'>
    y>|<cell|y<rprime|'> >|<cell| x>|<cell|
    y>|<cell|1>>>>><matrix|<tformat|<table|<row|<cell|f<rsub|11>>|<cell|f<rsub|12>>|<cell|f<rsub|13>>|<cell|f<rsub|21>>|<cell|f<rsub|22>>|<cell|f<rsub|23>>|<cell|f<rsub|31>>|<cell|f<rsub|32>>|<cell|
    f<rsub|33>>>>>><rsup|T>=0
  </equation*>

  each row can also be written as the fatten version of
  <math|x<rprime|'><rsup|T>x>.

  the solution shold be non-zero and up to scale, so 8 points are needed,
  solve can be obtained by SVD, use the sigular vector with the minimum
  sigular value as the solution.\ 

  By epipolar constrain, determinant of fundamental matrix is always 0 i.e.
  singular. So after the solution is found, we SVD the fundamental matrix \ 

  <\equation*>
    F<rprime|'>=U S V<rsup|T>
  </equation*>

  at last, reverse the normalization process, obtain the final result.

  <\equation*>
    F=N<rsup|T><rsub|1> F<rprime|''> N<rsub|2>
  </equation*>

  <subsubsection|Projection Matrix Estimation (DLT)>

  six coorespondence yields the projection matrix(3*4).

  get projection matrix

  <\equation*>
    x=<matrix|<tformat|<table|<row|<cell|x>>|<row|<cell|y>>|<row|<cell|1>>>>>=<matrix|<tformat|<twith|table-min-cols|4>|<table|<row|<cell|p<rsub|11>>|<cell|p<rsub|12>>|<cell|p<rsub|13>>|<cell|p<rsub|14>>>|<row|<cell|p<rsub|21>>|<cell|p<rsub|22>>|<cell|p<rsub|23>>|<cell|p<rsub|24>>>|<row|<cell|p<rsub|31>>|<cell|p<rsub|32>>|<cell|p<rsub|33>>|<cell|p<rsub|34>>>>>><matrix|<tformat|<table|<row|<cell|X>>|<row|<cell|Y>>|<row|<cell|Z>>|<row|<cell|1>>>>>=P
    X
  </equation*>

  <\equation*>
    0=x H x=x H P X=<wide|x|^>P X
  </equation*>

  <\equation*>
    <big|sum><rsup|i=3 j=4><rsub|i=0 j=0>a<rsub|i j> p<rsub|i j> =0\ 
  </equation*>

  <\equation*>
    <matrix|<tformat|<twith|table-min-cols|4>|<table|<row|<cell|a<rsub|11>>|<cell|a<rsub|12>>|<cell|a<rsub|13>>|<cell|a<rsub|14>>>|<row|<cell|a<rsub|21>>|<cell|a<rsub|22>>|<cell|a<rsub|23>>|<cell|a<rsub|24>>>|<row|<cell|a<rsub|31>>|<cell|a<rsub|32>>|<cell|a<rsub|33>>|<cell|a<rsub|34>>>>>>=<wide|x|^>X<rsup|T>
  </equation*>

  <subsubsection|Triangulation (DLT)>

  Recover the coordinate of a point in the world coordiante system from two
  calibrated camera and th projection coordinates of the point on each image
  plane. The constrains can be written as:\ 

  <\equation*>
    <choice|<tformat|<table|<row|<cell|x<rsub|l>=P<rsub|l>
    X<rsub|l>>>|<row|<cell|x<rsub|r>=P<rsub|r> X<rsub|r>>>>>>
  </equation*>

  A linear equation can be obtained.

  <\equation*>
    <matrix|<tformat|<table|<row|<cell|>>>>> P=A P=0
  </equation*>

  and the solution of P can be found by SVD the matrix <math|A<rsup|T>A>.

  <subsubsection|Standard Stereo Rectification>

  Rectification is the process of changing the relative pose of a pair of
  cameras into standard stereo configuration. It is essentially a homography
  of the original images. The three steps are: First send the epipoles to
  infinity by rotating both cameras by <math|R<rsub|rect>>.\ 

  <\equation*>
    R<rsub|rect>=<around*|[|e<rsub|1> e<rsub|2> e<rsub|3>|]><rsup|T>
  </equation*>

  Each row of <math|R<rsub|rect>> is obtained by the following algorithm.

  <\equation*>
    <choice|<tformat|<table|<row|<cell|e<rsub|1>=<dfrac|T|<around*|\<\|\|\>|T|\<\|\|\>>>>>|<row|<cell|e<rsub|2>=<dfrac|1|<sqrt|T<rsup|2><rsub|x>+T<rsup|2><rsub|y>>><around*|[|-T<rsub|y>
    T<rsub|x> 0|]><rsup|T>>>|<row|<cell|e<rsub|3>=e<rsub|1>\<times\>e<rsub|2>>>>>>
  </equation*>

  Second, rotate the reference camera by the relative roation of the other
  camera R so the epipolar lines became parallel. Finally, adjust the focal
  length so two image plane coincides.

  <subsection|RANSAC>

  RANSAC(RANdom SAmple Consensus) is a data filter used cross the entire
  system. In a set of samples with an underlying constrain, but contains both
  inliers and outliers, RANSAC attempts to select a small inlier-only subset
  by random selection within a finite number of trials and estimate this
  underlying constrain to further filter the dataset. RANSAC contains three
  major components:

  <\itemize>
    <item><code*|constrain> -- the underlying constrian of the samples

    <item><code*|constrain generator> -- the algorithm to estimate the
    constrain from a subset of samples

    <item><code*|dataset> -- the ``dirty'' data set
  </itemize>

  And three parameters:

  <\itemize>
    <item><code*|tolerance> -- tolerated error

    <item><code*|threshold><nbsp>-- estimated outlier percentage

    <item><code*|N>-- quota of trials
  </itemize>

  The process of RANSAC is a finite loop. Until no trials left, it randomly
  choose a subset of samples to estimate the underlying constrain e.g.
  fundamental matrix in two-view matching, then test the entire sample set,
  and calculate the percentage of rejected samples, if its greater than the
  threshold, this randomly selected subset may contain outlier which causes
  the estimated constrain to be inaccurate, and retry. If its less than the
  threshold, this estimated constrain is considered to be accuate, so we
  break out the loop and filter the dataset with this estimated constrain,
  and filter out the outliers using the estimated constrain.

  <subsection|Bundle Adjustment>

  Levenberg-Marqurdt Algorithm (LM) is the strategy we choose to perform
  non-linear optimization. Its the process to find the optimal solution so
  the difference between result and the measurement vector is minimized:

  <\equation*>
    <math-it|minimize> S<around*|(|p|)>=<around*|\<\|\|\>|f<around*|(|p|)>-y|\<\|\|\>>=<around*|\<\|\|\>|\<b-epsilon\><rsub|p>|\<\|\|\>>
  </equation*>

  In our case, we always use bundle adjustment to minimize error, so the
  measurement vector is assumed to be zero.

  In a small neighborhood, delta of the measurement is approximate to the
  prodect of Jacobian matrix and the step vector <math|\<b-delta\><rsub|p>> :

  <\equation*>
    f<around*|(|p+\<b-delta\><rsub|p>|)>\<thickapprox\>f<around*|(|p|)>+J\<b-delta\><rsub|p>
  </equation*>

  from which we obtian the equation:

  <\equation*>
    J<rsup|T>J\<b-delta\><rsub|p>=J<rsup|T>\<b-epsilon\><rsub|>
  </equation*>

  where

  <\equation*>
    N=J<rsup|T>J\<approx\>Hessian<around*|(|<around*|\<\|\|\>|\<b-epsilon\><rsub|p>|\<\|\|\>>|)>
  </equation*>

  What make the LMA special is it uses damping, the damped equation:

  <\equation*>
    <around*|(|N +\<b-mu\>I|)>\<b-delta\><rsub|p>=J<rsup|T>\<b-epsilon\><rsub|>
  </equation*>

  when far from optiamal, the damping value is large and the algorithm
  bahaves like steepest gradient descent, and when close to the optimal
  solution, damping value became smaller and behaves like Guass-Newton
  method. The damping value is updated during the iteration based on the
  current solution. The inital paramters are:

  <\equation*>
    \<b-mu\><rsub|0>=max<around*|(|N<rsub|ii>|)>\<nocomma\> , \<b-lambda\>=2
  </equation*>

  For each iteration, if the solution with damp value
  <math|\<b-mu\>/\<b-lambda\>> improved the solution, then the current
  optimal is updated and <math|\<b-mu\>/\<b-lambda\>> became the default damp
  value. Otherwise, if the current default damp value yields better solution,
  then the current optimal is updated and the damp value is unchaged. When
  both <math|\<b-mu\>/\<b-lambda\>> and <math|\<b-mu\>> do not yield better
  solution, than multiply the current damp value by <math|\<b-lambda\>> until
  the equation yields better solution than the current optimal.

  <section|Extract Image Features>

  There are several desirable feature extractors, we choose SIFT because the
  features it generates are scale/orientation invariant, and most existing
  SfM implementation uses SIFT as features extractor and yield disirable
  results.\ 

  <subsection|Scale Invariant>

  Features obtained from SIFT is scale invariant, which means features can be
  found and matched regardless what scale it is in. To accomplish that, we
  need to find feature points in LoGs (Laplacian of Guassian) of the image,
  due to the complexity of computing the second derivative of an image, we
  use DoGs(Difference of Guassian) to approximate. This is acomplished by
  using a series of octaves and scales of the original image, they simulates
  the image in various scales. We use 4 octave times 5 scales as Lowe's
  <cite|sift> paper suggested.

  To construct octave space, we directly shrink size of the image by half,
  each pixel is obtained by averaging the grayscale value of its window.In
  each octave, scalespace are constructed by convolution of the image with a
  series of a guassian kernels. The sigma of the guassian kernel is based on
  the relation between sigma of LoG and DoG.

  <\eqnarray*>
    <tformat|<table|<row|<cell|L<rsub|\<b-sigma\>>=\<nabla\><rsup|2>G<rsub|\<b-sigma\>>=G<rsub|\<b-sigma\>1>-G<rsub|\<b-sigma\>2>>|<cell|>|<cell|\<b-sigma\><rsub|1>=\<b-sigma\>/<sqrt|2>
    , \<b-sigma\><rsub|2>=<sqrt|2>\<b-sigma\>>>>>
  </eqnarray*>

  Sigma of guassian kernels in each octave is propotional by the factor of 2
  to satisfy the euqation above. After scalespaces are obtained, we can
  approximate LoG using DoG (Difference of Guassian) between adjacent scales.

  Afterwards, several filters are applied to choose high quality feature
  points. Contract filter \ will only accept maxium and minimum points in the
  DoG pyramid i.e. <math|3\<times\>3\<times\>3> window around the point.
  Corner filter will only accept corner by comparing the derivative of the
  point on each direction. Finally, scale invariant feature points are
  selected.\ 

  <subsection|Orientation Invariant>

  After finding scale inavriant feature points, we assign each feature point
  with an orientation from its neighborhood, then assign a discriptor to the
  point according to its orientation and neighborhood, so the descriptor is
  self-oriented, can be matched regardless the orientation of \ images. All
  following process are done in the DoG image where the point is found.\ 

  To calculate the orientation, we collect the gradients inside a window
  around the feature point and construct a histogram of gradient
  orientations, where orientations are divided into 36 discrete bins. All
  directions over 80% of the maxium is considered an valid orientation, which
  means one feaure point might coorespond to several oriented feature points.

  Finally, to calculate the descriptor, we collect gradient informations from
  a <math|16\<times\>16> window around the point aligned with its
  orientation. Gradient dirction are guass-weighted and relative to the
  orientation of feaure point. The window is then divided into a
  <math|4\<times\>4> grid of <math|4\<times\>4> sub-windows. In each
  sub-window, collect the histogram of gradient orientations with 8 bins.
  Afterwards, it will yield <math|16\<times\>8=128> numbers as the
  descripter. Finally, achive illumination invariant by normalizing the
  descriptor.

  <subsection|Usages>

  SIFT will generate a huge amount of feature points, it will be used in both
  camera registration phase and stereopsis phase. But only a small subset of
  features will be used in the camera registration phase, it's because of
  camera registratration need features that are globally unique, so it can be
  used to esitimate camera pose and parameters.\ 

  <section|Camera Registration>

  Camera registration i.e. camera calibration is the process of estimating
  camera parameters. Calibration matrix K, rotation matrix R and translation
  vector t are the three parameters we need to obtain in the calibration
  process.

  Intrinsic calibration is the process of obtaining the calibraton matrix K.
  The primary variable of K is focal length in pixel. Ideally we can obtain
  focal length in pixel from EXIF information, but in practice, focal length
  are given by mm, and most EXIF tag does not have CCD size or alternative
  attributes which allow us to calculate focal length in pixel, so we have to
  esitimate it by track coorespondence.

  Extrinsic calibration will recover the rotation matrix R and translation
  vector t. They can be extracted from both projection matrix and essential
  matrix, depend on the choice of algorithm.

  Camera parameters are very sensitive in view geometry, it will affect the
  estimated cooridante of every point it observes, and the relative pose
  bwtween cameras, consequently affect the entire system. Although bundle
  adjustment can refine the parameters, but sparse bundle adjustment only
  promises local optimization, the inital solution is still critical.\ 

  On the other hand, each image yields up to 10,000 features, but camera
  calibration algorithms require only a very small amount of matches e.g.
  DLT, 5-point-algorithm, so outliers in the input may affetct the result
  dramatically.

  So, from all points the feature extractor provides, choose matches with
  high percision, high contrast and globally unique is a crutial part of
  camera registration, and those points are called tracks, A track represents
  a distinct point in the scene, it have following constrains:

  <\itemize-dot>
    <item>Satisfy all geometric constrains

    <item>Does not have multiply matches in a single view i.e distinct
  </itemize-dot>

  <\algorithm*>
    <samp|inliers> := input tracks

    <samp|cameras> := input cameras

    <strong|do>:\ 

    <\indent>
      spsrse bundle adjustment (<samp|cameras>, <samp|inliers>)

      <strong|for each> camera <samp|vi>:

      <\indent>
        <math|d<rsub|80>> := 80 percentile of reprojection error

        <samp|threshold> := <math|clamp<around*|(|d<rsub|80>,4,16|)>>

        <samp|<math|outliers <rsub|vi> > >:= tracks with reprojection error
        larger then <samp|threshold>
      </indent>

      <samp|outliers> := <samp|<math|<big|cup>outliers<rsub| vi>>>

      <samp|inliers> := <samp|inliers> - <samp|outliers>
    </indent>

    <strong|while> <samp|outliers> is not empty
  </algorithm*>

  <subsection|Match Tracks>

  To select the desirable tracks, first match features for each pair of
  images. Instead using a threshold of the euclidean distance of two
  descriptor to indicate a match, we use a nearest neighbor strategy. When
  comparing two inages, for each feature in image 1, we find the 1st and 2nd
  nearest feature in image 2 and compare the ratio of the euclidean distance,
  \ if less then the ratio threshold, then accept as a match, which means
  it's a match and the only match in this image. then, use view constrains to
  filter the matches. for each image pair, use two view constrain,
  fundamental matrx. use RANSAC to estimate the fundamental matrix and filter
  outlier matches.

  To satisfy the second requirment, we construct image connectivity graph
  from all feature matches, tracks that contains multiple matches on a single
  image will be discarded. In this process, not only inconsistencies are
  filtered, we also merged two-view matches into a global image connectivity
  graph and global tracks, which can be used in next phase.

  <subsection|Register the First Pair of Cameras>

  The first pair of camera is crutial, we need to calibrate the inital pair
  of cameras and triangulate the inital set of tracks, it will be used in the
  incremental recovery phase to calibrate new cameras.

  We first guess the focal length in pixel. Then use 8-point-algorithm to
  estimate the fundamental matrix of the inital pair, and obtain its
  essential matrix using the guessed focal length, set one of the cmaera's
  local coordinate system as the world coordinate system i.e. set its
  rotation as I and translation as 0. Then carry out the other camera's
  extrinsic parameters from essential matrix. Then we select the tracks that
  are visible in this pair of cameras, and triangulate their world
  coordinates. Finally a global bundle adjustment is applied to refine the
  parameters, especially the focal length.

  <subsection|Incremetal Recovery>

  After the intial pair of cameras and tracks are registered, we can
  calibrate the rest of cameras one by one.\ 

  First find the next camera by selecting the one observes the most
  established tracks. Estimate its camera matrix from the coorespondence
  between each tracks's world coordinate and image coodinate using
  6-point-algorithm. Then a bundle adjustment is applied, only the new camera
  and tracks it observes can change. With camera matrix, we can obtain the
  calibration matrix, rotation matrix and translation vector using QR
  decomposition. Afterwards, triangulate tracks that are visible in the new
  camera. After each new track is added, a global bundle adjustment is
  applied. Repeat the process above until all cameras are registered or no
  camera observes enough tracks left.

  <section|Stereopsis>

  Stereopsis phase genterates a dense point cloud. It first locally
  reconstruct each two-view stereo, then merge them together with global
  constrians. Such as depth, visibility and epipolar constrain.

  <subsection|Normalized Cross Correlation (NCC)>

  <subsection|Multi-View Stereo>

  merge the two view stereos together. In two-view stereo matching process,
  points are triangulated within the two-view domain. From the matching lists
  we find all matches for a specific point and choose the optimal baseline
  then triangulate it again.

  <section|Application Framework>

  <subsection|Parallelization>

  There are two long-live threads:

  <\itemize-dot>
    <item>DOM Thread -- The default thread created by the browser, have full
    access to the DOM, the application logic of WebSFM is inside it.

    <item>WebSFM Thread -- SfM thread is a singleton thread created by the
    DOM thread, represents the serial logic of the SFM system.
  </itemize-dot>

  And three short-term threads

  <\itemize-dot>
    <item>SIFT feature generator thread

    <item>Two-View feature matching thread

    <item>Two-View stereo matching thread
  </itemize-dot>

  Multiple short term threads are created dynamically by the WebSFM thread
  according to the number of cores assigned to the application. They are
  implemented to be asynchronous. Each represents a worker for a task with no
  dynamic data dependency in the SFM system. For simplicity and performance
  concern, we only extracted three tasks to be asynchronous and parallel, all
  of them are well comfined in a specific scope -- a two-view pair or a
  single image.

  <subsection|Browser Related Details>

  There are some specific features that are crutial to the implementation of
  WebSFM, they provided us with a much more advanced programming environment
  in comparasion to the traditional client side javascript.\ 

  <\itemize-dot>
    <item>Typed Arrays -- Efficient data structure

    <item>2D Canvas -- \ Image Processing

    <item>IndexedDB -- Efficient data persistent

    <item>WebWorkers -- Parallel computation\ 

    <item>WebGL -- 3D Visualization
  </itemize-dot>

  <subsection|Samples>

  \;

  <small-figure||Two view matching using features generated by Lowe's SIFT>

  \;

  <small-figure||Dense point cloud generated by PMVS>

  <\bibliography|bib|tm-bibtex|thesis.bib>
    <\bib-list|4>
      <bibitem*|1><label|bib-sift>Distinctive image features from
      scale-invariant keypoints.<newblock> <with|font-shape|italic|IJCV>, ,
      2004.<newblock>

      <bibitem*|2><label|bib-lma>The design and implementation of a generic
      sparse bundle adjustment software package based on the
      levenberg-marquardt algorithm.<newblock> <with|font-shape|italic|>, ,
      2004.<newblock>

      <bibitem*|3><label|bib-modeltheworld>Modeling the world from internet
      photo collections.<newblock> <with|font-shape|italic|>, ,
      2007.<newblock>

      <bibitem*|4><label|bib-handheld>Visual modeling with a hand-held
      camera.<newblock> <with|font-shape|italic|>. <newblock>
    </bib-list>
  </bibliography>
</body>

<\initial>
  <\collection>
    <associate|page-medium|paper>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|2>>
    <associate|auto-10|<tuple|3.1.5|4>>
    <associate|auto-11|<tuple|3.1.6|4>>
    <associate|auto-12|<tuple|3.2|5>>
    <associate|auto-13|<tuple|3.3|5>>
    <associate|auto-14|<tuple|4|6>>
    <associate|auto-15|<tuple|4.1|6>>
    <associate|auto-16|<tuple|4.2|6>>
    <associate|auto-17|<tuple|4.3|6>>
    <associate|auto-18|<tuple|5|7>>
    <associate|auto-19|<tuple|5.1|7>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-20|<tuple|5.2|7>>
    <associate|auto-21|<tuple|5.3|7>>
    <associate|auto-22|<tuple|6|8>>
    <associate|auto-23|<tuple|6.1|8>>
    <associate|auto-24|<tuple|6.2|8>>
    <associate|auto-25|<tuple|7|8>>
    <associate|auto-26|<tuple|7.1|8>>
    <associate|auto-27|<tuple|7.2|8>>
    <associate|auto-28|<tuple|7.3|9>>
    <associate|auto-29|<tuple|2|9>>
    <associate|auto-3|<tuple|1|2>>
    <associate|auto-30|<tuple|3|9>>
    <associate|auto-31|<tuple|3|9>>
    <associate|auto-32|<tuple|3|?>>
    <associate|auto-33|<tuple|10.4|?>>
    <associate|auto-34|<tuple|10.5|?>>
    <associate|auto-35|<tuple|10.5|?>>
    <associate|auto-36|<tuple|10.5|?>>
    <associate|auto-4|<tuple|3|3>>
    <associate|auto-5|<tuple|3.1|3>>
    <associate|auto-6|<tuple|3.1.1|3>>
    <associate|auto-7|<tuple|3.1.2|3>>
    <associate|auto-8|<tuple|3.1.3|3>>
    <associate|auto-9|<tuple|3.1.4|4>>
    <associate|bib-handheld|<tuple|4|9>>
    <associate|bib-lma|<tuple|2|9>>
    <associate|bib-modeltheworld|<tuple|3|9>>
    <associate|bib-pmvspaper|<tuple|1|?>>
    <associate|bib-sift|<tuple|1|9>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      handheld

      modeltheworld

      lma

      sift

      sift
    </associate>
    <\associate|figure>
      <tuple|normal|Framework of WebSFM|<pageref|auto-3>>

      <tuple|normal|Two view matching using features generated by Lowe's
      SIFT|<pageref|auto-29>>

      <tuple|normal|Dense point cloud generated by PMVS|<pageref|auto-30>>
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Introduction>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Overview>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Mathematical
      Context> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>View Geometry
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|2tab>|3.1.1<space|2spc>Coordinate System
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|2tab>|3.1.2<space|2spc>View Constrains
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7>>

      <with|par-left|<quote|2tab>|3.1.3<space|2spc>Normalized Eight Point
      Algorithm (DLT) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8>>

      <with|par-left|<quote|2tab>|3.1.4<space|2spc>Projection Matrix
      Estimation (DLT) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9>>

      <with|par-left|<quote|2tab>|3.1.5<space|2spc>Triangulation (DLT)
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|2tab>|3.1.6<space|2spc>Standard Stereo
      Rectification <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>RANSAC
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Bundle Adjustment
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|4<space|2spc>Extract
      Image Features> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14><vspace|0.5fn>

      <with|par-left|<quote|1tab>|4.1<space|2spc>Scale Invariant
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|1tab>|4.2<space|2spc>Orientation Invariant
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16>>

      <with|par-left|<quote|1tab>|4.3<space|2spc>Usages
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|5<space|2spc>Camera
      Registration> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.5fn>

      <with|par-left|<quote|1tab>|5.1<space|2spc>Match Tracks
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19>>

      <with|par-left|<quote|1tab>|5.2<space|2spc>Register the First Pair of
      Cameras <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20>>

      <with|par-left|<quote|1tab>|5.3<space|2spc>Incremetal Recovery
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|6<space|2spc>Stereopsis>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22><vspace|0.5fn>

      <with|par-left|<quote|1tab>|6.1<space|2spc>Normalized Cross Correlation
      (NCC) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <with|par-left|<quote|1tab>|6.2<space|2spc>Multi-View Stereo
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|7<space|2spc>Application
      Framework> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25><vspace|0.5fn>

      <with|par-left|<quote|1tab>|7.1<space|2spc>Parallelization
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <with|par-left|<quote|1tab>|7.2<space|2spc>Browser Related Details
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27>>

      <with|par-left|<quote|1tab>|7.3<space|2spc>Samples
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-31><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>