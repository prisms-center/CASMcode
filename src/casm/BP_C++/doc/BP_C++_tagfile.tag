<?xml version='1.0' encoding='ISO-8859-1' standalone='yes' ?>
<tagfile>
  <compound kind="file">
    <name>BP_Comb.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___comb_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <class kind="class">BP::BP_Comb</class>
    <namespace>BP</namespace>
  </compound>
  <compound kind="file">
    <name>BP_coords.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p__coords_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <includes id="_b_p___parse_8h" name="BP_Parse.h" local="yes" imported="no">BP_Parse.h</includes>
    <includes id="_b_p__z_parse_8h" name="BP_zParse.h" local="yes" imported="no">BP_zParse.h</includes>
    <class kind="class">BP::ucs_coord_class</class>
    <class kind="class">BP::frac_coord_class</class>
    <class kind="class">BP::cart_coord_class</class>
    <namespace>BP</namespace>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>adf83378bf2df54eb11f3e3d73a17bbcb</anchor>
      <arglist>(string &amp;s, ucs_coord_class &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate_x</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6883df46588b374f3a2d9d9c25047564</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, double theta)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate_y</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a1cce1ddbce46e7d44f6d12f85ac43263</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, double theta)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate_z</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a907e903cd0a21a9330538f17710986c2</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, double theta)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a61473d2bd8d303e2799880def1e4db0d</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, BP_Vec&lt; BP_Vec&lt; double &gt; &gt; &amp;R)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>frac_to_cart</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a8eb2891968d2814be56345dd98fef5f1</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;sc_v, frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>cart_to_frac</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac846f4551ef5763af0e49bca7394b489</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;sc_v, cart_coord_class c)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_Geo.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___geo_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <includes id="_b_p___g_vec_8h" name="BP_GVec.h" local="yes" imported="no">BP_GVec.h</includes>
    <includes id="_b_p__useful_8h" name="BP_useful.h" local="yes" imported="no">BP_useful.h</includes>
    <class kind="class">BP::V_Element</class>
    <class kind="class">BP::CH_data_class</class>
    <class kind="class">BP::DT_data_class</class>
    <class kind="class">BP::Geo</class>
    <namespace>BP</namespace>
    <member kind="function">
      <type>bool</type>
      <name>vector_is</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a8bf92cfec78b6521dae37bee14cc109c</anchor>
      <arglist>(const VectorXd &amp;v, double val, double tol)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>vector_is</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>adbbf19be50c05cad04919f6383f9c379</anchor>
      <arglist>(const VectorXd &amp;v, const VectorXd &amp;v2, double tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_bins</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7473ce465a007c06933049c9d4a94837</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Gen_GVec_Member * &gt; &gt; &amp;point_bins, BP_GVec&lt; V_Element &gt; &amp;pts, int curr_facet)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_GVec.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___g_vec_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <class kind="class">BP::BP_Gen_GVec</class>
    <class kind="class">BP::BP_GVec</class>
    <class kind="class">BP::BP_Group_Data_class</class>
    <class kind="class">BP::BP_Graph_Data_class</class>
    <class kind="class">BP::BP_Gen_GVec_Member</class>
    <class kind="class">BP::BP_GVec_Member</class>
    <class kind="class">BP::BP_Gen_Group</class>
    <class kind="class">BP::BP_Group</class>
    <class kind="class">BP::BP_Edge_Data_class</class>
    <class kind="class">BP::BP_Gen_Vertex</class>
    <class kind="class">BP::BP_Gen_Edge</class>
    <class kind="class">BP::BP_Gen_Graph</class>
    <class kind="class">BP::BP_Graph</class>
    <namespace>BP</namespace>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a518124e0130b45970250906b6c3cedcf</anchor>
      <arglist>(BP_GVec&lt; T &gt; &amp;gvec, const T &amp;i1)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_Parse.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___parse_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <class kind="class">BP::BP_Parse</class>
    <class kind="class">BP::BP_bParse</class>
    <class kind="class">BP::BP_Write</class>
    <class kind="class">BP::BP_bWrite</class>
    <namespace>BP</namespace>
    <member kind="function">
      <type>void</type>
      <name>cut_start_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a0799283ad09a64045a7957a0a08525a0</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_end_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2634a00015e7af65c76d1d8455433081</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_start_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aa6b12b4a9e7061b28c1c24c1eb5efbff</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_end_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a756904ad1701941183b583d15ecbc6bb</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>itos</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad5570112aa5bb6bf9b28444f72ba9d08</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2d1f18ac9249f10ad9e3396387a25d28</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a80626dbfaeb80c4c6a663f7687d959dc</anchor>
      <arglist>(double i1, int prec)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos_fixed</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad5fc760a57b5723c438d1f168640cfdf</anchor>
      <arglist>(double i1, int prec)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos_scientific</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aea238465f9d5e28402d382cb359d5276</anchor>
      <arglist>(double i1, int prec)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>stoi</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a5a5bbcf75843ba43794b9e5eee1d6503</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>stod</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a23fa80f4b0f2010abaa6fadb40ccc6ce</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>stoi</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a18766df52bedbcd8942d05e704423062</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; double &gt;</type>
      <name>stod</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6dc80213cb070819905d1900bd963bb1</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>any_match</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a86189f854dae939133df07220a9bec5c</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;str_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a0529803bf780d3f524d7adebd7d4ae48</anchor>
      <arglist>(string &amp;s1, const string &amp;com)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a054ca39b107c320a5e21a7507a13a62e</anchor>
      <arglist>(string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ab5db361afaa9926c0320d99771cc2e71</anchor>
      <arglist>(const string &amp;s1, const string &amp;com)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac496467838cdc0e14517cd82c658d60b</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>trim</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a40adc5871069fbdcfab1102ea3ea4093</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;delim_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aa720faab9e121579f5daada357717911</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a55bde3f6b088f5ffa196c19e0bbea78f</anchor>
      <arglist>(string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a57c335fd1b7460d765d7f0bccd872b7a</anchor>
      <arglist>(BP_Vec&lt; string &gt; &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a1183cc0ee959ccdbd8de9f0cd39c2d9b</anchor>
      <arglist>(BP_Vec&lt; string &gt; &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a37666591102c644345cf32c2b4fbae4f</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7095aecfffb1d0c70544b2aae0fe4ba9</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aa59cbaa1876bfd2b7afe30608eb4b23a</anchor>
      <arglist>(const string &amp;s1, const string &amp;delim)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af0e8c2a26bc1b49d0a6007af3653ebf1</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;delim_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a5dbdc5b7f521efb94b4893bf7cbaebe4</anchor>
      <arglist>(const string &amp;s1, const string &amp;delim)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aae5f169e49794a4780103f01ff36b3d4</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;delim_list)</arglist>
    </member>
    <member kind="function">
      <type>long int</type>
      <name>stoli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>abd30106d2f5cce341bc4651dd74cf0fb</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>stouli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a9b540242526f2624e1873ea863529cac</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; long int &gt;</type>
      <name>stoli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aef501a56d74d657aa342c5e5f40bbf28</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; unsigned long int &gt;</type>
      <name>stouli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a9ad0f610a8c66fadfda48ad70dc3ca3f</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>next_string</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a8254ab32acbae23984134547459ce15a</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>next_double</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7f755ed5145444f599f8b13f0b745532</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>next_int</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac264a4d64ff74177b780054798e67e0a</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>next_bool</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a4e6712211345ab18e80209ad3d159620</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a3736c5c3a5a9f0a9379dcb54b4e1585a</anchor>
      <arglist>(string &amp;s, int &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af818eb11458ac055e08e17d39361b8e9</anchor>
      <arglist>(string &amp;s, double &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ace7639f5214a3bc8fbf902c72e3c9630</anchor>
      <arglist>(string &amp;s, string &amp;i)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_Plot.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___plot_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <includes id="_b_p___parse_8h" name="BP_Parse.h" local="yes" imported="no">BP_Parse.h</includes>
    <class kind="class">BP::BP_RGB</class>
    <class kind="class">BP::BP_Text</class>
    <class kind="class">BP::BP_Plot_Data</class>
    <class kind="class">BP::BP_Axes</class>
    <class kind="class">BP::BP_Plot</class>
    <namespace>BP</namespace>
  </compound>
  <compound kind="file">
    <name>BP_RVG.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___r_v_g_8h</filename>
    <includes id="_b_p___g_vec_8h" name="BP_GVec.h" local="yes" imported="no">BP_GVec.h</includes>
    <includes id="_b_p__useful_8h" name="BP_useful.h" local="yes" imported="no">BP_useful.h</includes>
    <class kind="class">BP::BP_RVG_base</class>
    <class kind="class">BP::BP_RVG_linear</class>
    <class kind="class">BP::BP_RVG_tree</class>
    <namespace>BP</namespace>
  </compound>
  <compound kind="file">
    <name>BP_StopWatch.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___stop_watch_8h</filename>
    <class kind="class">BP::BP_StopWatch</class>
    <namespace>BP</namespace>
  </compound>
  <compound kind="file">
    <name>BP_unix_cheat_sheet.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p__unix__cheat__sheet_8h</filename>
    <member kind="function">
      <type>echo test mail s Hello foo bar com saving without saving changes i insert esc exit whatever mode you are in move to beginning of line move to end of line nG move to nth line of the file G move to last line of the file H move to top of screen M move to middle of screen L move to bottom of screen move to</type>
      <name>associated</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>aec7ed957095140cabd88cbc26e457ef8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>svnadmin create path to repository svn import project</type>
      <name>file</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>ac6685c268746b7a302f6397ef79b133b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>svnadmin create path to repository svn import project</type>
      <name>change2</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>ad89d0ec0966690ec48ad5b8d652f38a8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>svnadmin create path to repository svn import project etc svn update bash script sh file1 file2 echo echo ls grep CAR</type>
      <name>VAR</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>a9afca942140457795830f136af5d4ec9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>do ls</type>
      <name>$i</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>a32f1d3c2e2b29456c786fcf466430bee</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>echo test mail s Hello foo bar com</type>
      <name>commands</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>a1fbce94760e66314494228f5825c2cf6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>echo test mail s Hello foo bar com saving</type>
      <name>changes</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>ada643a79f641a74336c87a8873ff9846</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Ctrl u move up page Ctrl d move down page Ctrl b move up page Ctrl f move down page dd delete current line ndd delete n lines yy copy line p paste string search for string string search backwards for string n go to next instance of string N go to previous instance of</type>
      <name>string</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>a3bc965dfcce84081b92d73c96529d140</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>Ctrl u move up page Ctrl d move down page Ctrl b move up page Ctrl f move down page dd delete current line ndd delete n lines yy copy line p paste string search for string string search backwards for string n go to next instance of string N go to previous instance of Ctrl M standard</type>
      <name>variables</name>
      <anchorfile>_b_p__unix__cheat__sheet_8h.html</anchorfile>
      <anchor>a36c502a27600a738aff3a4d4c2625a37</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_useful.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p__useful_8h</filename>
    <namespace>BP</namespace>
    <member kind="function">
      <type>void</type>
      <name>BP_pause</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7ba7983014d0d66b7dbeffaa3a3d7d12</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sqr</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aac3721fb09390af0f08b276ee142bbb6</anchor>
      <arglist>(double a)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_even</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a820c94bf5dc5c77f537d68e71fc78e79</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>int_pow</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6eb6af11b0923730accf4150b714137c</anchor>
      <arglist>(int x, int p)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>ulint_pow</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>acb85bc49f04876e0912bc1382554949c</anchor>
      <arglist>(unsigned long int x, unsigned long int p)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_Vasp.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___vasp_8h</filename>
    <includes id="_b_p___parse_8h" name="BP_Parse.h" local="yes" imported="no">BP_Parse.h</includes>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <includes id="_b_p__coords_8h" name="BP_coords.h" local="yes" imported="no">BP_coords.h</includes>
    <class kind="class">BP::BP_POSCAR_class</class>
    <namespace>BP</namespace>
  </compound>
  <compound kind="file">
    <name>BP_Vec.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p___vec_8h</filename>
    <includes id="_b_p__useful_8h" name="BP_useful.h" local="yes" imported="no">BP_useful.h</includes>
    <class kind="class">BP::BP_Vec</class>
    <namespace>BP</namespace>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad7ee0e02647ce29cab85abf64a774c33</anchor>
      <arglist>(ostream &amp;sout, const BP_Vec&lt; U &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sum</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aab23fabd49fc2a059109727d620a0099</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a4bf76b6269ad65ea48d3d5c611026568</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af7d70a5ddce69ce89d2c8fa15cb6bbff</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, unsigned long int &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7c3fdecca80146909c1c93f38d95c0c5</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6db95421fa6009f51b864fc50142ae72</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aba2e68a95c4e1534c119549ea0de236f</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, unsigned long int &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a81b54d847dc5b43aaea8dd48a4b6ec18</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sum</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aee1775fe60bb1a8fc4cf7ea9665e9503</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7d56e68b1e0979264f11641f922e0f8b</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a15528f2bb0c9dfb8509cfafe7414cda7</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, unsigned long int &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7d71228af45cd7f07c114fe6057bbb86</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6b9541745404662ecb71a8a3c4ba50ef</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2eb2529b05b23ba354798180b42096b4</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, unsigned long int &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a0b0f0db079addd98046600f83aeadc7e</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7975db3f02124d8db7f2dffc30269228</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rms</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad0aae4e69af388eb04b8daa36955e8eb</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mag</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ae63bd02ec8497193c522b61fbb8deb98</anchor>
      <arglist>(const BP_Vec&lt; T &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ab098028012ad690827fd2085459c6655</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const T &amp;i1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af7cb5a0cbe9775651c0059440efb243d</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const T &amp;i1, double eps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac6248634a4b28ea8b1cd4cef2688d7d3</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const BP_Vec&lt; T &gt; &amp;vec2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>mem_size</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a163ccd82ce29a1fa7773b0bd1b8e7f14</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>mem_size</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a41ca84a56730f092fe04d3d701951f0d</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>mem_size</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af06f76feee18494d4e8f3343def178c9</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a234026e463eef2f4a7c81483f398e75e</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a75b6e575680db93ed1f9fbeafbc86da7</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a792607dfd145ac826ca640af80f1f1e5</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2b2ac3461bc8de5db377f350df8d9879</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a198ffc0841305a3eb9cd680ebfeff257</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>BP_zParse.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_b_p__z_parse_8h</filename>
    <includes id="_b_p___vec_8h" name="BP_Vec.h" local="yes" imported="no">BP_Vec.h</includes>
    <includes id="_b_p___parse_8h" name="BP_Parse.h" local="yes" imported="no">BP_Parse.h</includes>
    <class kind="class">BP::BP_zParse</class>
    <class kind="class">BP::BP_bzParse</class>
    <class kind="class">BP::BP_zWrite</class>
    <class kind="class">BP::BP_bzWrite</class>
    <namespace>BP</namespace>
  </compound>
  <compound kind="file">
    <name>MersenneTwister.h</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_mersenne_twister_8h</filename>
    <class kind="class">MTRand</class>
    <member kind="function">
      <type>std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>_mersenne_twister_8h.html</anchorfile>
      <anchor>a059061d50a1e54ee3067d4e1dbdd7c64</anchor>
      <arglist>(std::ostream &amp;os, const MTRand &amp;mtrand)</arglist>
    </member>
    <member kind="function">
      <type>std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>_mersenne_twister_8h.html</anchorfile>
      <anchor>a45b02a702835a3be42171c5c2dc79b2d</anchor>
      <arglist>(std::istream &amp;is, MTRand &amp;mtrand)</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>MTexample.cpp</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>_m_texample_8cpp</filename>
    <includes id="_mersenne_twister_8h" name="MersenneTwister.h" local="yes" imported="no">MersenneTwister.h</includes>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>_m_texample_8cpp.html</anchorfile>
      <anchor>af3ed9c200de85b53c94cd18764b246a2</anchor>
      <arglist>(int argc, char *const argv[])</arglist>
    </member>
  </compound>
  <compound kind="file">
    <name>template.cpp</name>
    <path>/Users/cmg/BP_C++/</path>
    <filename>template_8cpp</filename>
    <member kind="function">
      <type>int</type>
      <name>main</name>
      <anchorfile>template_8cpp.html</anchorfile>
      <anchor>ae66f6b31b5ad750f1fe042a706a4e3d4</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>MTRand</name>
    <filename>class_m_t_rand.html</filename>
    <member kind="enumvalue">
      <name>N</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>ab8fea37d16b55e1a0fe06149e325f1b6a60f472facea8fabd42765cd91273db7b</anchor>
      <arglist></arglist>
    </member>
    <member kind="enumvalue">
      <name>SAVE</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a7d9f4f1783a4e45f7834dd5174dfc2a1a3899803ea0d4da3018d311ed4902d9cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="typedef">
      <type>unsigned long</type>
      <name>uint32</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a45478edf9e24dcd2a5164bac3889d6a2</anchor>
      <arglist></arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MTRand</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a172bc7e7cf1e578ef3f9c90a8cee3eb1</anchor>
      <arglist>(const uint32 oneSeed)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MTRand</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a380e79e0192b46426abcefa6e2dd082e</anchor>
      <arglist>(uint32 *const bigSeed, uint32 const seedLength=N)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MTRand</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a265dc65546e26073c0d5f8787b045a1d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>MTRand</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>aff69d4a4ec88475bab03a295e8fb0f60</anchor>
      <arglist>(const MTRand &amp;o)</arglist>
    </member>
    <member kind="function">
      <type>uint32</type>
      <name>randInt</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>ad1008efd4fe0e8aae30459c2c58cfe35</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>uint32</type>
      <name>randInt</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a3515bbf6e1b46680a4ce6968451942b6</anchor>
      <arglist>(const uint32 n)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rand</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a76d129a2d850c24ff4a0613f299cf3a5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rand</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>aa4fe82fc27fd81414ce7554093a9766b</anchor>
      <arglist>(const double n)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>randExc</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>afd05e468983b3a3d66ce0f403bd666af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>randExc</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>aa1e89d6c7ac8737567b3ccf8fe70b6de</anchor>
      <arglist>(const double n)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>randDblExc</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a4d3a475aa72fe6d1a6d7d9e16d6a732e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>randDblExc</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a1a81d8f00de8f553d4b8626d64e1c544</anchor>
      <arglist>(const double n)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>operator()</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>abbb87a08d622d58fdee0eea4cb5471a0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rand53</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a15f4daf79febbe4ff43c3e6ce2c4fcbe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>randNorm</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a4c284f626b6d40a0367ff2a949ea1944</anchor>
      <arglist>(const double mean=0.0, const double stddev=1.0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>seed</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a1e21a79e0a30225fffe924229e34a923</anchor>
      <arglist>(const uint32 oneSeed)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>seed</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a5758103776b131e8ea46b6dc1b9fb267</anchor>
      <arglist>(uint32 *const bigSeed, const uint32 seedLength=N)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>seed</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>ad88ea3363d55bafb62826bbd130279c2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>save</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>ad60e0f3f5c90baab75b74f9a2ccae871</anchor>
      <arglist>(uint32 *saveArray) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>load</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a8302e9a8cd16d8dfc536a85bf2f68be0</anchor>
      <arglist>(uint32 *const loadArray)</arglist>
    </member>
    <member kind="function">
      <type>MTRand &amp;</type>
      <name>operator=</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a3a6eb21add6f6ef4ce2d3280f2518521</anchor>
      <arglist>(const MTRand &amp;o)</arglist>
    </member>
    <member kind="enumvalue">
      <name>M</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a10c3437be98225f5b0beee1ed8c033c8a133070000b798889cd75535ea0d5bb71</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>initialize</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a9b9a20998f5c805af6301ce5c37dcfc3</anchor>
      <arglist>(const uint32 oneSeed)</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>void</type>
      <name>reload</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a1d5fcb69d83f4d2fd653883c8352f86c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="protected">
      <type>uint32</type>
      <name>hiBit</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a45eea926a0602e4bb5c0b90b04779826</anchor>
      <arglist>(const uint32 u) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>uint32</type>
      <name>loBit</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a6f5a4a532e1c3acd42052046594205be</anchor>
      <arglist>(const uint32 u) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>uint32</type>
      <name>loBits</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>ad846f81f7abfc1b20c51d1563b8e5d45</anchor>
      <arglist>(const uint32 u) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>uint32</type>
      <name>mixBits</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>abdd5587252ed1ac89cb274e4bf4881da</anchor>
      <arglist>(const uint32 u, const uint32 v) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>uint32</type>
      <name>magic</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a8539a48116c85704c5101981cb0823e7</anchor>
      <arglist>(const uint32 u) const </arglist>
    </member>
    <member kind="function" protection="protected">
      <type>uint32</type>
      <name>twist</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>acf32530212717166e3d02dd3cc0b68c4</anchor>
      <arglist>(const uint32 m, const uint32 s0, const uint32 s1) const </arglist>
    </member>
    <member kind="function" protection="protected" static="yes">
      <type>static uint32</type>
      <name>hash</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a486885d03f38c844315d002e6312fa23</anchor>
      <arglist>(time_t t, clock_t c)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>uint32</type>
      <name>state</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a2c87f537429bf0b0f6a452c22b9eebba</anchor>
      <arglist>[N]</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>uint32 *</type>
      <name>pNext</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a2b80858137c88fe69d4d2bdc665bcf93</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>left</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a98eabf568c88f121e44f487397f32495</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend std::ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a059061d50a1e54ee3067d4e1dbdd7c64</anchor>
      <arglist>(std::ostream &amp;os, const MTRand &amp;mtrand)</arglist>
    </member>
    <member kind="friend">
      <type>friend std::istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_m_t_rand.html</anchorfile>
      <anchor>a45b02a702835a3be42171c5c2dc79b2d</anchor>
      <arglist>(std::istream &amp;is, MTRand &amp;mtrand)</arglist>
    </member>
  </compound>
  <compound kind="namespace">
    <name>BP</name>
    <filename>namespace_b_p.html</filename>
    <class kind="class">BP::BP_Comb</class>
    <class kind="class">BP::ucs_coord_class</class>
    <class kind="class">BP::frac_coord_class</class>
    <class kind="class">BP::cart_coord_class</class>
    <class kind="class">BP::V_Element</class>
    <class kind="class">BP::CH_data_class</class>
    <class kind="class">BP::DT_data_class</class>
    <class kind="class">BP::Geo</class>
    <class kind="class">BP::BP_Gen_GVec</class>
    <class kind="class">BP::BP_GVec</class>
    <class kind="class">BP::BP_Group_Data_class</class>
    <class kind="class">BP::BP_Graph_Data_class</class>
    <class kind="class">BP::BP_Gen_GVec_Member</class>
    <class kind="class">BP::BP_GVec_Member</class>
    <class kind="class">BP::BP_Gen_Group</class>
    <class kind="class">BP::BP_Group</class>
    <class kind="class">BP::BP_Edge_Data_class</class>
    <class kind="class">BP::BP_Gen_Vertex</class>
    <class kind="class">BP::BP_Gen_Edge</class>
    <class kind="class">BP::BP_Gen_Graph</class>
    <class kind="class">BP::BP_Graph</class>
    <class kind="class">BP::BP_Parse</class>
    <class kind="class">BP::BP_bParse</class>
    <class kind="class">BP::BP_Write</class>
    <class kind="class">BP::BP_bWrite</class>
    <class kind="class">BP::BP_RGB</class>
    <class kind="class">BP::BP_Text</class>
    <class kind="class">BP::BP_Plot_Data</class>
    <class kind="class">BP::BP_Axes</class>
    <class kind="class">BP::BP_Plot</class>
    <class kind="class">BP::BP_RVG_base</class>
    <class kind="class">BP::BP_RVG_linear</class>
    <class kind="class">BP::BP_RVG_tree</class>
    <class kind="class">BP::BP_StopWatch</class>
    <class kind="class">BP::BP_POSCAR_class</class>
    <class kind="class">BP::BP_Vec</class>
    <class kind="class">BP::BP_zParse</class>
    <class kind="class">BP::BP_bzParse</class>
    <class kind="class">BP::BP_zWrite</class>
    <class kind="class">BP::BP_bzWrite</class>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>adf83378bf2df54eb11f3e3d73a17bbcb</anchor>
      <arglist>(string &amp;s, ucs_coord_class &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate_x</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6883df46588b374f3a2d9d9c25047564</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, double theta)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate_y</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a1cce1ddbce46e7d44f6d12f85ac43263</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, double theta)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate_z</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a907e903cd0a21a9330538f17710986c2</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, double theta)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt; &amp;</type>
      <name>rotate</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a61473d2bd8d303e2799880def1e4db0d</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;v_list, BP_Vec&lt; BP_Vec&lt; double &gt; &gt; &amp;R)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>frac_to_cart</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a8eb2891968d2814be56345dd98fef5f1</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;sc_v, frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>cart_to_frac</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac846f4551ef5763af0e49bca7394b489</anchor>
      <arglist>(BP_Vec&lt; cart_coord_class &gt; &amp;sc_v, cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>vector_is</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a8bf92cfec78b6521dae37bee14cc109c</anchor>
      <arglist>(const VectorXd &amp;v, double val, double tol)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>vector_is</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>adbbf19be50c05cad04919f6383f9c379</anchor>
      <arglist>(const VectorXd &amp;v, const VectorXd &amp;v2, double tol)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print_bins</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7473ce465a007c06933049c9d4a94837</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Gen_GVec_Member * &gt; &gt; &amp;point_bins, BP_GVec&lt; V_Element &gt; &amp;pts, int curr_facet)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a518124e0130b45970250906b6c3cedcf</anchor>
      <arglist>(BP_GVec&lt; T &gt; &amp;gvec, const T &amp;i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_start_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a0799283ad09a64045a7957a0a08525a0</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_end_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2634a00015e7af65c76d1d8455433081</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_start_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aa6b12b4a9e7061b28c1c24c1eb5efbff</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_end_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a756904ad1701941183b583d15ecbc6bb</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>itos</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad5570112aa5bb6bf9b28444f72ba9d08</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2d1f18ac9249f10ad9e3396387a25d28</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a80626dbfaeb80c4c6a663f7687d959dc</anchor>
      <arglist>(double i1, int prec)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos_fixed</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad5fc760a57b5723c438d1f168640cfdf</anchor>
      <arglist>(double i1, int prec)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>dtos_scientific</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aea238465f9d5e28402d382cb359d5276</anchor>
      <arglist>(double i1, int prec)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>stoi</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a5a5bbcf75843ba43794b9e5eee1d6503</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>stod</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a23fa80f4b0f2010abaa6fadb40ccc6ce</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>stoi</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a18766df52bedbcd8942d05e704423062</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; double &gt;</type>
      <name>stod</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6dc80213cb070819905d1900bd963bb1</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>any_match</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a86189f854dae939133df07220a9bec5c</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;str_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a0529803bf780d3f524d7adebd7d4ae48</anchor>
      <arglist>(string &amp;s1, const string &amp;com)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a054ca39b107c320a5e21a7507a13a62e</anchor>
      <arglist>(string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ab5db361afaa9926c0320d99771cc2e71</anchor>
      <arglist>(const string &amp;s1, const string &amp;com)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>cut_comment</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac496467838cdc0e14517cd82c658d60b</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>trim</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a40adc5871069fbdcfab1102ea3ea4093</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;delim_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aa720faab9e121579f5daada357717911</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a55bde3f6b088f5ffa196c19e0bbea78f</anchor>
      <arglist>(string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a57c335fd1b7460d765d7f0bccd872b7a</anchor>
      <arglist>(BP_Vec&lt; string &gt; &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a1183cc0ee959ccdbd8de9f0cd39c2d9b</anchor>
      <arglist>(BP_Vec&lt; string &gt; &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>clean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a37666591102c644345cf32c2b4fbae4f</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;com_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7095aecfffb1d0c70544b2aae0fe4ba9</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aa59cbaa1876bfd2b7afe30608eb4b23a</anchor>
      <arglist>(const string &amp;s1, const string &amp;delim)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize_whtspace</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af0e8c2a26bc1b49d0a6007af3653ebf1</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;delim_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a5dbdc5b7f521efb94b4893bf7cbaebe4</anchor>
      <arglist>(const string &amp;s1, const string &amp;delim)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>tokenize</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aae5f169e49794a4780103f01ff36b3d4</anchor>
      <arglist>(const string &amp;s1, const BP_Vec&lt; string &gt; &amp;delim_list)</arglist>
    </member>
    <member kind="function">
      <type>long int</type>
      <name>stoli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>abd30106d2f5cce341bc4651dd74cf0fb</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>stouli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a9b540242526f2624e1873ea863529cac</anchor>
      <arglist>(const string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; long int &gt;</type>
      <name>stoli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aef501a56d74d657aa342c5e5f40bbf28</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; unsigned long int &gt;</type>
      <name>stouli</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a9ad0f610a8c66fadfda48ad70dc3ca3f</anchor>
      <arglist>(const BP_Vec&lt; string &gt; &amp;s_list)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>next_string</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a8254ab32acbae23984134547459ce15a</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>next_double</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7f755ed5145444f599f8b13f0b745532</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>next_int</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac264a4d64ff74177b780054798e67e0a</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>next_bool</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a4e6712211345ab18e80209ad3d159620</anchor>
      <arglist>(string &amp;s1)</arglist>
    </member>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a3736c5c3a5a9f0a9379dcb54b4e1585a</anchor>
      <arglist>(string &amp;s, int &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af818eb11458ac055e08e17d39361b8e9</anchor>
      <arglist>(string &amp;s, double &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>string &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ace7639f5214a3bc8fbf902c72e3c9630</anchor>
      <arglist>(string &amp;s, string &amp;i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>BP_pause</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7ba7983014d0d66b7dbeffaa3a3d7d12</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sqr</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aac3721fb09390af0f08b276ee142bbb6</anchor>
      <arglist>(double a)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_even</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a820c94bf5dc5c77f537d68e71fc78e79</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>int_pow</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6eb6af11b0923730accf4150b714137c</anchor>
      <arglist>(int x, int p)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>ulint_pow</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>acb85bc49f04876e0912bc1382554949c</anchor>
      <arglist>(unsigned long int x, unsigned long int p)</arglist>
    </member>
    <member kind="function">
      <type>ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad7ee0e02647ce29cab85abf64a774c33</anchor>
      <arglist>(ostream &amp;sout, const BP_Vec&lt; U &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>sum</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aab23fabd49fc2a059109727d620a0099</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a4bf76b6269ad65ea48d3d5c611026568</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af7d70a5ddce69ce89d2c8fa15cb6bbff</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, unsigned long int &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7c3fdecca80146909c1c93f38d95c0c5</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6db95421fa6009f51b864fc50142ae72</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aba2e68a95c4e1534c119549ea0de236f</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, unsigned long int &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a81b54d847dc5b43aaea8dd48a4b6ec18</anchor>
      <arglist>(const BP_Vec&lt; int &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>sum</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>aee1775fe60bb1a8fc4cf7ea9665e9503</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7d56e68b1e0979264f11641f922e0f8b</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a15528f2bb0c9dfb8509cfafe7414cda7</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, unsigned long int &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7d71228af45cd7f07c114fe6057bbb86</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;min_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a6b9541745404662ecb71a8a3c4ba50ef</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2eb2529b05b23ba354798180b42096b4</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, unsigned long int &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a0b0f0db079addd98046600f83aeadc7e</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list, BP_Vec&lt; unsigned long int &gt; &amp;max_index)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mean</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a7975db3f02124d8db7f2dffc30269228</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>rms</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ad0aae4e69af388eb04b8daa36955e8eb</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;i_list)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mag</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ae63bd02ec8497193c522b61fbb8deb98</anchor>
      <arglist>(const BP_Vec&lt; T &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ab098028012ad690827fd2085459c6655</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const T &amp;i1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af7cb5a0cbe9775651c0059440efb243d</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const T &amp;i1, double eps)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_once</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>ac6248634a4b28ea8b1cd4cef2688d7d3</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const BP_Vec&lt; T &gt; &amp;vec2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>mem_size</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a163ccd82ce29a1fa7773b0bd1b8e7f14</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>mem_size</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a41ca84a56730f092fe04d3d701951f0d</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>mem_size</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>af06f76feee18494d4e8f3343def178c9</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a234026e463eef2f4a7c81483f398e75e</anchor>
      <arglist>(BP_Vec&lt; T &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a75b6e575680db93ed1f9fbeafbc86da7</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a792607dfd145ac826ca640af80f1f1e5</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a2b2ac3461bc8de5db377f350df8d9879</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>namespace_b_p.html</anchorfile>
      <anchor>a198ffc0841305a3eb9cd680ebfeff257</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; BP_Vec&lt; T &gt; &gt; &gt; &gt; &gt; &amp;vec, const T &amp;val)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Comb</name>
    <filename>class_b_p_1_1_b_p___comb.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_Comb</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>ab1ab7753b28d79e7455cf9f1871a9c73</anchor>
      <arglist>(int n, int k)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a1a17d4898a9b5464da44b80d6ad8b3f3</anchor>
      <arglist>(int n, int k)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>total_combs</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a632229721f8d90b9faa44cbe0745af3e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a3913eddbd54bbebb78f8f7cf217db447</anchor>
      <arglist>(int i) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>acfbf8968b08fd06cb3f354ec6906395f</anchor>
      <arglist>(int i) const </arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>get_bit_string</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a7c7beab10f3c873a940904a13edf8aff</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; bool &gt;</type>
      <name>get_all</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a6a239784d9ae617be0e79fa86278e8dd</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>increment</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>ada4ea5d85ffe5a7c9fdd6be940c9aa68</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_N</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a8c89f215fe98942cefe12d0db8b0efee</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_K</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a8b2ba2539e6df0093e036e31a5925a6e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>aa06bfda21d42fac8557edb07a6efd3a7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_count</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a1810d321375dbdf1328c417467291336</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>complete</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>acc650a373a9379a2b09d8d0faef3e09b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>fix</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a4b5820c7a1d780ae4cd2b3bdc12eb75f</anchor>
      <arglist>(int i1, bool i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unfix</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>afe062ffe2e06518a47d80e543167c7a3</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>unfix_all</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a8e6494aca82537f35d157c418998f8ad</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; bool &gt;</type>
      <name>val_list</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>aa0aef9b8135318c348430b2d59ae0672</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; unsigned int &gt;</type>
      <name>pos_list</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>aae6917c4faf1599239c25e1a680d3ac2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; unsigned int &gt;</type>
      <name>fix_pos_list</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a1ffd50a4112e2964cbad3eadaa280489</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; bool &gt;</type>
      <name>fix_val_list</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>aee571a21685ddac05e4e1eb4b3f33d8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>curr_index</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a40c962ac755416f251802f0d9028a59b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>is_complete</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>adfdea5ee1bc2f7249beb8ffb93ef389a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>N</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a0fed3b5f1ae2dcef1c98dfd5afa78741</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>K</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a6dd4e03d9859994cab6d2ab3c476cfba</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned long int</type>
      <name>count</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>a6de681ab4ff2bac968a8379d9d2e1f28</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1_b_p___comb.html</anchorfile>
      <anchor>af5372b2196ad6706c275b921ccadd65f</anchor>
      <arglist>(ostream &amp;outstream, const BP_Comb &amp;c)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::ucs_coord_class</name>
    <filename>class_b_p_1_1ucs__coord__class.html</filename>
    <member kind="function">
      <type></type>
      <name>ucs_coord_class</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a4d305f9b7da8d0493b2b00faa96bde75</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>ucs_coord_class</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>aec59f98f8848111dccb9a9a392c425dc</anchor>
      <arglist>(int i0, int i1, int i2, int i3)</arglist>
    </member>
    <member kind="function">
      <type>int &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>ad32e228cafd1cdb8ebc2ce9655b6fd39</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>const int &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>aaff2148aeba131362cd9c4e52090324c</anchor>
      <arglist>(int i1) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator==</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>af3a6dd90ed0df1506beae5bda9a15ef0</anchor>
      <arglist>(ucs_coord_class B) const </arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>operator!=</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a52bf35bb0a149b406361769af3ee55cb</anchor>
      <arglist>(ucs_coord_class B) const </arglist>
    </member>
    <member kind="function">
      <type>ucs_coord_class</type>
      <name>operator+</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a65322972e528a638b1f6a4a468b15fd3</anchor>
      <arglist>(ucs_coord_class B)</arglist>
    </member>
    <member kind="function">
      <type>ucs_coord_class</type>
      <name>operator+</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a54938789159372c3e39cb9b8bc9a51de</anchor>
      <arglist>(ucs_coord_class B) const </arglist>
    </member>
    <member kind="function">
      <type>ucs_coord_class</type>
      <name>operator-</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>ad148fc9850a860a45fcab579fa8e5480</anchor>
      <arglist>(ucs_coord_class B)</arglist>
    </member>
    <member kind="function">
      <type>ucs_coord_class</type>
      <name>operator-</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a5fa4cddb838396e5e60605bc33cb9053</anchor>
      <arglist>(ucs_coord_class B) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a10e1a27b7932be60934aa291ff9fa558</anchor>
      <arglist>(BP_bWrite &amp;fout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>read</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>aa74a7bcd7d9e5bb9ac61ab12be4795b9</anchor>
      <arglist>(BP_bParse &amp;fin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a66a8180d9ee0e25631ec7d1fcd9c49b5</anchor>
      <arglist>(BP_bzWrite &amp;fout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>read</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a6d2978bb9f7aff7e5a410e94cfd8835d</anchor>
      <arglist>(BP_bzParse &amp;fin)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a7bc6f15a19c2cdce31767dbf915694a2</anchor>
      <arglist>(int i0, int i1, int i2, int i3)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>coord</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a7785a6dac742e0d8ea312a844bb19bae</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a613486a09e296d07da7c8c60a242d668</anchor>
      <arglist>(ostream &amp;outstream, const ucs_coord_class &amp;c)</arglist>
    </member>
    <member kind="friend">
      <type>friend istream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_b_p_1_1ucs__coord__class.html</anchorfile>
      <anchor>a4f921b7a68efc9ad8c431d8c8cf4c17f</anchor>
      <arglist>(istream &amp;instream, ucs_coord_class &amp;c)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::frac_coord_class</name>
    <filename>class_b_p_1_1frac__coord__class.html</filename>
    <member kind="function">
      <type></type>
      <name>frac_coord_class</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a977209bd2c091ffd765903ab80ae4d20</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>frac_coord_class</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a8041fe37292b13923f63514ba32acc69</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>double &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a0cdb9a3f566b2764c122a407e54f9488</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>ad833385e4f2225821ca8606f2ccd38c8</anchor>
      <arglist>(int i1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a30663c601fc6068a1691eda4b77c1fd3</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>operator+</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a6a135f65ede863ce7d98559229ac528d</anchor>
      <arglist>(frac_coord_class f) const </arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>operator-</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>adcfaff91f7a5a2a6e923fbd420457a22</anchor>
      <arglist>(frac_coord_class f) const </arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>af32198ead0c1cb1b3eb0fcd3e43580d0</anchor>
      <arglist>(double d) const </arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>operator/</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>af5b07561e7b365f20fbe19f3ba8ab5b9</anchor>
      <arglist>(double d) const </arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class &amp;</type>
      <name>operator+=</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>ac09c246ba8cefe379fd7e11e57d86576</anchor>
      <arglist>(frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class &amp;</type>
      <name>operator-=</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a2716197b5d949e3728442a419e6c5d73</anchor>
      <arglist>(frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>equal</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>aad4c1a3113bf0915533feec8c562e127</anchor>
      <arglist>(frac_coord_class f, double eps)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mag</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a78cbc0d61bc0bdaa1c05a9aa8ea881b0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>in_unit</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>a442cc64f14a2d6cf9e866c437cede3f3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>frac_coord_class</type>
      <name>in_unit</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>aeac01de408a6a99506745ded16307cd3</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>coord</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>aaa179383b00d056c9f19d5facb22cc43</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1frac__coord__class.html</anchorfile>
      <anchor>ae2f155cf9d256cba668dbdfdbaa0096b</anchor>
      <arglist>(ostream &amp;outstream, const frac_coord_class &amp;c)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::cart_coord_class</name>
    <filename>class_b_p_1_1cart__coord__class.html</filename>
    <member kind="function">
      <type></type>
      <name>cart_coord_class</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>aef04dbe3e3d5ea7d8b5df22360160d09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>cart_coord_class</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a921509b9f150ac9cea431e12f8f498b2</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>double &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>aa60ca2bddfd36f788a98d6b45c0f5148</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a18f1956d1351005b1db8c8adabc9c1d9</anchor>
      <arglist>(int i1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>abd07aadef4cd05e7160e837fea18bdba</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>mag</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>ad7a070b90213c0af3094cfee345def05</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>norm</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a128a1fb29d894364a13e0e4712f149c5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>dot</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a20cb8eb18eda6cfd24fd645d8b64d86b</anchor>
      <arglist>(const cart_coord_class &amp;c)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>cross</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a7a54c0b7defa832d6736dd99c4c21cb2</anchor>
      <arglist>(const cart_coord_class &amp;b)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator+</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>ac225db91b1caf5d280a40004a0d44424</anchor>
      <arglist>(cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator-</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a88a731184f834c34c3a2605161561f2f</anchor>
      <arglist>(cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class &amp;</type>
      <name>operator+=</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a4ca3251ef38a19f25c51abf1f944df77</anchor>
      <arglist>(cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class &amp;</type>
      <name>operator-=</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a6ef9adf1c519a59dfac1dfcab4198ad0</anchor>
      <arglist>(cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a5ae107ff08c492893912f63ebf5cd95e</anchor>
      <arglist>(int x)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a0ca3c49b406d6bc4672d4b9bdc9957c6</anchor>
      <arglist>(int x) const </arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a37556a2407bb21cc12fa3c56e384cb2c</anchor>
      <arglist>(unsigned long int x)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a2385516136444f93402a8e09256a23c6</anchor>
      <arglist>(unsigned long int x) const </arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a97e68fdec2903330b4f9dcb7d6e43c48</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator*</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a7a0c516a7e4e3d07b0153ea3945bd5b4</anchor>
      <arglist>(double x) const </arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator/</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a4ac2f4fc323a5dfcbddf811e40d2fa29</anchor>
      <arglist>(double x)</arglist>
    </member>
    <member kind="function">
      <type>cart_coord_class</type>
      <name>operator/</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a45fa28501195aa24ded689ada13cea31</anchor>
      <arglist>(double x) const </arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>coord</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a6c513ca9488be3aaf3fdd609d774037f</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1cart__coord__class.html</anchorfile>
      <anchor>a3e3d3a514328064349c93b8634f38666</anchor>
      <arglist>(ostream &amp;outstream, const cart_coord_class &amp;c)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::V_Element</name>
    <filename>class_b_p_1_1_v___element.html</filename>
    <member kind="function">
      <type>V_Element</type>
      <name>intersect</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a52fc0604d59016c93ab168854c34fe4f</anchor>
      <arglist>(const V_Element &amp;E1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>V_Element</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a9c34591cbe52ff6c30c4611fdfae165c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a20f59d39cc0b5e07cddad5077333b1a1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>V_Element</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a426f185b53fad17b2e90d324f2d9044b</anchor>
      <arglist>(int i1, VectorXd &amp;i2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>V_Element</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a416ea649ab1f797346969a0dd4a3e6b0</anchor>
      <arglist>(int i1, VectorXd &amp;i2, VectorXd &amp;i3)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>V_Element</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>aafbca37e4ee027104b6b152b4f39f35a</anchor>
      <arglist>(int i1, VectorXd &amp;i2, VectorXd &amp;i3, VectorXd &amp;i4)</arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>dim</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a4c747d8602c9e24b8092590ef169f381</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>VectorXd</type>
      <name>pos</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a8a268b50ae063896eb147286203a121f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>rank</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>ad507a0d7f0ffb81b1171188fb4b2b8da</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>VectorXd</type>
      <name>vertex_mean</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>af25f5a0e6c66cf6094300aa1e3557b6f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>VectorXd</type>
      <name>out_norm</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a1c44adebb2e8c037d6c33fa8c07be1df</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>nvolume</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>af48d79faa6fe43588e83a233ec3d53e0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>cut</name>
      <anchorfile>class_b_p_1_1_v___element.html</anchorfile>
      <anchor>a3d0a5744326ec23438dc9be1920e0c37</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::CH_data_class</name>
    <filename>class_b_p_1_1_c_h__data__class.html</filename>
    <member kind="function">
      <type></type>
      <name>CH_data_class</name>
      <anchorfile>class_b_p_1_1_c_h__data__class.html</anchorfile>
      <anchor>ab5d755b7d87263474645493429cb9c89</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>CH_data_class</name>
      <anchorfile>class_b_p_1_1_c_h__data__class.html</anchorfile>
      <anchor>a76e415eb73ea30694a603157b1f007e7</anchor>
      <arglist>(bool i1, double i2, const VectorXd &amp;i3)</arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>is_hull_point</name>
      <anchorfile>class_b_p_1_1_c_h__data__class.html</anchorfile>
      <anchor>a49484eacf058f68db2f722215caa19d8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Vec&lt; int &gt;</type>
      <name>closest_facet</name>
      <anchorfile>class_b_p_1_1_c_h__data__class.html</anchorfile>
      <anchor>a847a6c621f4a3e419334c48439f30ab3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>double</type>
      <name>dist_to_hull</name>
      <anchorfile>class_b_p_1_1_c_h__data__class.html</anchorfile>
      <anchor>a5b3c0e3d9defd3b077a90cc42a5fce5b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>VectorXd</type>
      <name>vec_to_hull</name>
      <anchorfile>class_b_p_1_1_c_h__data__class.html</anchorfile>
      <anchor>aef77fd38a13752e658e11cc67e0e6c5a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::DT_data_class</name>
    <filename>class_b_p_1_1_d_t__data__class.html</filename>
    <member kind="function">
      <type></type>
      <name>DT_data_class</name>
      <anchorfile>class_b_p_1_1_d_t__data__class.html</anchorfile>
      <anchor>a6e00a188ea1531b252ae92286d9f7275</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>DT_data_class</name>
      <anchorfile>class_b_p_1_1_d_t__data__class.html</anchorfile>
      <anchor>a54db29e4cc6d782a4560b18636e89c06</anchor>
      <arglist>(const VectorXd &amp;i1)</arglist>
    </member>
    <member kind="variable">
      <type>VectorXd</type>
      <name>ccenter</name>
      <anchorfile>class_b_p_1_1_d_t__data__class.html</anchorfile>
      <anchor>a29888fc7b533768acefd446c3cc512da</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::Geo</name>
    <filename>class_b_p_1_1_geo.html</filename>
    <member kind="function">
      <type></type>
      <name>Geo</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a16e9b15c70508d9a1776075db37c2dcb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Geo</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a14c333adeb0f9e29c1378dd74b045873</anchor>
      <arglist>(const Geo &amp;geo)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>Geo</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a03165c6b4adee85c6f066522509d31d2</anchor>
      <arglist>(const MatrixXd &amp;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~Geo</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ac0ffc246048a48855ac1c338766806af</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a4dc7619e08d9f1b320b4a80516b99326</anchor>
      <arglist>(const MatrixXd &amp;)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>get_dim</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a845fa69aa5ba9687b9abfac889cfed21</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ae797bb94fd3ba4b2e67aefda8555550e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>pos</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a20eedaaa33eae85eeba9e2d904e80554</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_verbosity</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aaeb58e11ef1cd22dbe0e732449153532</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; BP_Vec&lt; int &gt; &gt;</type>
      <name>get_equivalent_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a40f9a97dc21c13badc3fadfbc8d32e25</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_equivalent_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a55bcda1149fc634dde8cc73f7cb3cf00</anchor>
      <arglist>(ostream &amp;sout)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_Geo_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a3337ef4e2f96a16d80dfffa12632842a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_CH_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8937af6a4990c406955b776f4d425260</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_VD_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad3dbe38afdbb7da388cd3f3691b85669</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_Geo_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a26dd748c1a9f18a4bd42c5c584d7f260</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_CH_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a3b2542e6402d01ff72c2eb00008ae2b3</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_VD_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a34d7b5536230a298bed96f29a2841f75</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>calc_CH</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ac514bedcf068e7fa6063d9e234367d73</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CH_bottom</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a606846494d260f59e512d6eb9522fdfe</anchor>
      <arglist>(const VectorXd &amp;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CH_wholehull</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a0c909d0b58d53450038132888bcd2d82</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>CH_volume</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a30374103cd013b23a9343f8c4b3db78d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>CH_area</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a1cc56702b3a0b9ec7920bd12d91dc3aa</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>CH_write</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa7cc2ad1c714007c8684b09da8ece2b9</anchor>
      <arglist>(ostream &amp;)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>CH_closest_facet</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a9dd69dc621025459833aa8074cedffe6</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>CH_dist_to_hull</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa80d2ab37f63458d6ed88ef0115bfc1d</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>CH_vec_to_hull</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>abc5cd4b8db5a73e516cf7a3c2e151f84</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>CH_dist_to_hull</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a541ded6d1326d201eb04132e7e9e2a0e</anchor>
      <arglist>(const VectorXd &amp;)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>CH_facets_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>acab578542935e5f4fc936006c1fd5ed8</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>CH_facets_area</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a779c6fa4c1dc2dc02bdc0e3ad5e11adc</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>CH_facets_norm</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ada9ae5cd6fe2d0b5b28fe70725c2fd8d</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>CH_facets_nborverts_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a3938604c2ac9ec8de314e77954716d9b</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>CH_facets_nborverts</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa9ea20eb1facca74a04c450a7b8a0624</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>CH_facets_nborfacets_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a0b4b8dfc3dfbc667e642c678c1e90b81</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>CH_facets_nborfacets</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a4eaa30d6b2c477ee2acbd55a1e4b2c4c</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>CH_verts_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8e1cefa77af62b931139ec109f317a8e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>CH_verts_pos</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>af83d0501be6cf66043f75f2d6cad54a2</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>CH_verts_indices</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a2f659b8d20a3716f0e1a4105737aee7b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>CH_verts_nborverts_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a00813665d7fd7f282c73f346ff9f518a</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>CH_verts_nborverts</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad1edeee494c4f8e9d86577258e33f2fc</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>CH_verts_nborfacets_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aea1c0c31be5106dca360ae53273f21ef</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>CH_verts_nborfacets</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a9f41ee3a26accb946ef4a627b3e76ee9</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>calc_DT</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ab09a266e9c63498f501fe57214d23969</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>DT_write</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a9abff716a8b433ce48ed1f26ea5d9605</anchor>
      <arglist>(ostream &amp;)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>DT_volume</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a6a0a441e2569512df0fa3aaf563a7930</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DT_simps_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad6627cd96d881dd33ced9cc0c2f3e85b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>DT_simps_volume</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad917677ee33a90af9619b62cf2051cf2</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>DT_simps_circumcenter</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad5c2491259ec38b81a7dc384988f412c</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DT_simps_nborverts_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>abe0ec3d46c73f44aad15e2a4e6a5d812</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>DT_simps_nborverts</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ae316abe8a94ec3f81ce9b82987d32b38</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DT_simps_nborsimps_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a49c598bf8e7f13ae8be6694357686c5e</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>DT_simps_nborsimps</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aefba9755ff72ec5f2187398d83a028e2</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DT_verts_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8792b3eabac55cc46da0d7a92a6642a0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>DT_verts_pos</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ab917af8c8b418848a06ad9adc3c3cff2</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DT_verts_nborverts_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad7fd0849d5c69e1b69b4df3ced5e6304</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>DT_verts_nborverts</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a2895f4d0df47eb6b5909d24c51449e09</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>DT_verts_nborsimps_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a9a17a73435eda0cc77474b8e8b121f6f</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>DT_verts_nborsimps</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a85f09dea45f708b166316fa3f547d512</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>calc_VD</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a1eb3dd130791b961db91fc7acf69fb21</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>VD_elements_size</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a0990d356dea88707195a67d6167f7866</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>VD_elements_verts</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a0118de7280b22cfdaaecb9944c23c0a1</anchor>
      <arglist>(int, int)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>VD_elements_volume</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a454188e5dd05e01b2b1fc2c2339549c2</anchor>
      <arglist>(int, int)</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>VD_vertex_pos</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a2e7d5f7afa96208219c3a3ca3ac9a055</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>VectorXd</type>
      <name>VD_elements_vertex_mean</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa39ca91769c0fb86879d2543a94a0408</anchor>
      <arglist>(int, int)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>generate_CH</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa4b506c327481c31e290dfb8a02727e2</anchor>
      <arglist>(BP_GVec&lt; V_Element &gt; &amp;, BP_Group&lt; V_Element &gt; &amp;, BP_GVec&lt; V_Element &gt; &amp;, BP_Graph&lt; V_Element, bool &gt; &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>form_simplex_from_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ae7401831ea4f56ffa894cd4ac66b83b6</anchor>
      <arglist>(BP_GVec&lt; V_Element &gt; &amp;, BP_Group&lt; V_Element &gt; &amp;, BP_GVec&lt; V_Element &gt; &amp;, BP_Graph&lt; V_Element, bool &gt; &amp;, BP_Group&lt; V_Element &gt; &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>check_hull_point</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>adc5d7e7db51c833a31fde8f015fbed76</anchor>
      <arglist>(const VectorXd &amp;, BP_GVec&lt; V_Element &gt; &amp;, BP_Graph&lt; V_Element, bool &gt; &amp;, BP_Vec&lt; BP_Gen_Vertex * &gt; &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>check_hull_point_BFS</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a06ca18160170d93e99a2f01478327090</anchor>
      <arglist>(const VectorXd &amp;, BP_GVec&lt; V_Element &gt; &amp;, BP_Graph&lt; V_Element, bool &gt; &amp;, BP_Vec&lt; BP_Gen_Vertex * &gt; &amp;, int)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>BP_Gen_Vertex *</type>
      <name>form_facet_from_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8d6aa4a500cf977a1aedeee9c0757e33</anchor>
      <arglist>(BP_Gen_Vertex *, BP_Vec&lt; BP_Gen_Vertex * &gt; &amp;, BP_GVec&lt; V_Element &gt; &amp;, BP_Graph&lt; V_Element, bool &gt; &amp;, const VectorXd &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>write_hull</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a792cd84ae68e8ff9b868f084e0b05297</anchor>
      <arglist>(int, int, BP_GVec&lt; V_Element &gt; &amp;, BP_Group&lt; V_Element &gt; &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>bool</type>
      <name>use_facet</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ada34a6262aed11d20674b49c6746dcac</anchor>
      <arglist>(int)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>CH_calc_distances</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a69ae92fbce77769ef41da3696ae0c364</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>DT_calc_data</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a6e904d1e625850d1ca36783c5c1cbc5b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>VD_calc_data</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a69471726cedefbd2ea34f821e1a025f6</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>VD_face_finder</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a3eca6ec1a8f451c0ca0b595465e9631a</anchor>
      <arglist>(int, int, BP_Gen_GVec_Member *, BP_Vec&lt; int &gt; &amp;, BP_Vec&lt; int &gt; &amp;, BP_Vec&lt; int &gt; &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>BP_Gen_GVec_Member *</type>
      <name>VD_form_facet_from_ccenters</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ab415ff53ef1853cc0a3e74a628ca6f22</anchor>
      <arglist>(int, BP_GVec&lt; V_Element &gt; &amp;, BP_Vec&lt; int &gt; &amp;)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>int</type>
      <name>compare</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a485c4e083eb34abc44f38ba825dd5fc4</anchor>
      <arglist>(double, double, double)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>VectorXd</type>
      <name>scale</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a166d5b6ac2c89da9cec0c7e7bcb92122</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_GVec&lt; V_Element &gt;</type>
      <name>points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ae86ebd794fc2461beb59845adb4e8cdd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_GVec&lt; bool &gt;</type>
      <name>connections</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a2934b37fb37852c6c63ed7672812d956</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>dim</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a7e62c8a4167d52d1f4429cc1551cab2d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>use_bottom</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa405a4c9c25d290d76ca331186a26f0d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>VectorXd</type>
      <name>bottom</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>abba367656880e711618477bafcb6533f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>Geo_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ac80116ddccc190ef250b4dec47f83242</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>verbosity</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>adac49eb8bcb01674685b1ca7f2147b50</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_Vec&lt; int &gt; &gt;</type>
      <name>equivalent_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ab573ba383d476e71d5454d5b8ac53507</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>CH_exists</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa9f773843e4a7223416c107e7048da74</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>ch_area</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ad8fbaf1690aa2bd07d983d22ac82d1dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>ch_volume</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>abd4a2eacd214bc59a26f7f4193dcc9e1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>CH_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>acdb9674fbd7b93cbbf8c8a7d0e3a7f19</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Group&lt; V_Element &gt;</type>
      <name>CH_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a30e5fce0ccd56569f532f93a2fe3c142</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_GVec&lt; V_Element &gt;</type>
      <name>CH_facets</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a5fd3ee9a9240a9a927b655cada5ae8b4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Group&lt; V_Element &gt;</type>
      <name>CH_filter_points</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8469b70d8f386e7c5d24d8c7d85ba4ab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Group&lt; V_Element &gt;</type>
      <name>CH_filter_facets</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ab9376c60fb5ab2e3b32e5b7525382aa2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; CH_data_class &gt;</type>
      <name>CH_data</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a840c7e680ba53dab0f95ed4398b702d2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>CH_data_exists</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8dc5226d97595dd6939a88a564b83066</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>DT_exists</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a6f429303ac27f9d0d4487b1f21ee7234</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>dt_volume</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a7a879ea559b657927b08ea8025364958</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_GVec&lt; V_Element &gt;</type>
      <name>DT_simplexes</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a8b587025ec7489f87287c274c3fcf33e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; DT_data_class &gt;</type>
      <name>DT_data</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a49d1e0e9abffdec8c3cc2b90843f48cd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>DT_data_exists</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a9ce4b8d646adb2cceb4ca50288c79427</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>VD_exists</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>aa3c3376f2a89092696679aa752fa4e96</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>VD_tol</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>ac44d789f130ca74ce5f84b88de708c06</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_GVec&lt; V_Element &gt; *</type>
      <name>VD_elements</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>abc445fd606595c2952ef58b00cfaed88</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>VD_data_exists</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a602d3b5e99108cb357d857d107f9e2ee</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Graph&lt; V_Element, bool &gt;</type>
      <name>CH_graph</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>abff3f5385e0def968d3817a15be2c339</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Graph&lt; V_Element, bool &gt;</type>
      <name>DT_graph</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a021a78573819774aad13c4d334eeda99</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Graph&lt; V_Element, bool &gt;</type>
      <name>VD_graph</name>
      <anchorfile>class_b_p_1_1_geo.html</anchorfile>
      <anchor>a2061a5f32bfcb5465f6b16cd4a6c46a0</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Gen_GVec</name>
    <filename>class_b_p_1_1_b_p___gen___g_vec.html</filename>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>ad4d6d99a64359c296e3ff696d1d77b78</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>ae7c859aea30630fa0b6c035c982bc7aa</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ordered_remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a9700280f81e085ec8e671b439ffa1a2b</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a374c4c64c855c2180f74abf2c51de28c</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ordered_remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a2c30aefbe86601e7d25b12cbb5cec0fd</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>swap</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>ac6092d7aa8d8fb337cadd4a9cc4f4a74</anchor>
      <arglist>(unsigned long int, unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ordered_swap</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a1431e6ed0a6aff10bc4b0a53cd59abea</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a85d4f2ebd7c494398251e7c7431cd013</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a3a6f5e0a52dc4d768ce7b42d8b717ed3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>ac0749a41f641835bc4e2c87a700ed5ff</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>ab967570e6fbf14166fb2f8ce4abbe164</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a123fc8a88388746d9fd5832859196d96</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cap_increment</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a0e0d8b321c73928923dea46c7e1e4b7f</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>aa901c2f42804b98913c2c71727539e30</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Gen_GVec</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a75cff871c09d927db62b4f1e3a13b645</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_Gen_GVec</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a968fe20552ca78dd8c8ef9900d2a9811</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>N_max</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a33d5be0be94ca874cd1fde69e29be512</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>N_alloc</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>abe14822128fdbb2351921d3fd5fab4e5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>N</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>ad1cc943be5b23d74440908c94c79290d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>cap_incr</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>a8a11ad831feac859b056ae90441f8199</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>BP_Gen_GVec_Member **</type>
      <name>val</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec.html</anchorfile>
      <anchor>aa96ca412f1133f3372ce347b88534359</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_GVec</name>
    <filename>class_b_p_1_1_b_p___g_vec.html</filename>
    <templarg>T</templarg>
    <base>BP::BP_Gen_GVec</base>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>ae8734fbe48fab5f1bdde2a8a3ff858e9</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a181ec74eda052db0f4879fc3edaafde5</anchor>
      <arglist>(BP_Gen_Vertex *)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>acc921331d4c4c3a7f5d51704340caea1</anchor>
      <arglist>(BP_Gen_Edge *)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a64a92160c6cee89d0f811dace574628d</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>aa0530aa05d622bac45390a7bb8ded449</anchor>
      <arglist>(unsigned long int) const </arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a5bb63d2f9a9ddce42207c8ce62ea032b</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>aaaf3ce746ba37c18d850d0d799933b5d</anchor>
      <arglist>(BP_Gen_GVec_Member *) const </arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a946fad1004a30ef1983a5468f1ec9d12</anchor>
      <arglist>(BP_Gen_Vertex *)</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a422439a3cf3592cb679fc4d3e8da3291</anchor>
      <arglist>(BP_Gen_Vertex *) const </arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a3de61a9e284bfa050bbfaedb19e8c785</anchor>
      <arglist>(BP_Gen_Edge *)</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a69fddd82d53ec9284f95da617cea5a23</anchor>
      <arglist>(BP_Gen_Edge *) const </arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a46a24d73489d57cfdb6fedbb124c9f8a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a42ad108e704ca3f7489ab2c1955d0d52</anchor>
      <arglist>(const T &amp;)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a8415b5fab0f8d50a4fc75a053e3d2e5a</anchor>
      <arglist>(const U &amp;)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add_in_place</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a19ef1b2bd1b3c5329a5aec67e985767e</anchor>
      <arglist>(const U &amp;, unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>last</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>aa041aa3c41cac78481b8830535a9ae3c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>last</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>ac3abbe3d352b3eeb2844cdac54dd930d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec.html</anchorfile>
      <anchor>a5dbf85edaf90721cb108306b89a70a4c</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Group_Data_class</name>
    <filename>class_b_p_1_1_b_p___group___data__class.html</filename>
    <member kind="variable">
      <type>BP_Gen_Group *</type>
      <name>address</name>
      <anchorfile>class_b_p_1_1_b_p___group___data__class.html</anchorfile>
      <anchor>adc9f114878dd372720b0741997227ad8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned long int</type>
      <name>index</name>
      <anchorfile>class_b_p_1_1_b_p___group___data__class.html</anchorfile>
      <anchor>a3412b0b54a59234610d6e2effb6e3e35</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Graph_Data_class</name>
    <filename>class_b_p_1_1_b_p___graph___data__class.html</filename>
    <member kind="variable">
      <type>BP_Gen_Graph *</type>
      <name>address</name>
      <anchorfile>class_b_p_1_1_b_p___graph___data__class.html</anchorfile>
      <anchor>a6b87d5b33fa0fd6077d3f22bde4b171f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>int</type>
      <name>type</name>
      <anchorfile>class_b_p_1_1_b_p___graph___data__class.html</anchorfile>
      <anchor>a2fdcdccd2acb23efa2b4ac612cd2d4ab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned long int</type>
      <name>index</name>
      <anchorfile>class_b_p_1_1_b_p___graph___data__class.html</anchorfile>
      <anchor>ace786fe89067f7bfa34ba8c7f060f079</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Gen_GVec_Member</name>
    <filename>class_b_p_1_1_b_p___gen___g_vec___member.html</filename>
    <member kind="function">
      <type>BP_Gen_GVec *</type>
      <name>GVec_home</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>afd95414358b2454c44eb0ce8b9029406</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>GVec_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a887fa3ffaa32dc30f6910b74f01ea584</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_GVec_home</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a2c94b1a1a65fe03964dc2df3ad9e0bb7</anchor>
      <arglist>(BP_Gen_GVec *i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_GVec_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a9f643b32884370a48700ab5b881340e1</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>group_size</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>abc44d8494efc1d78af8f74b2b2b6c91d</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>graph_size</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>aef5e521c9d617ed6121ad19e1d92020e</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Group *</type>
      <name>group</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a7e72374b16e303490b4339e0a2481a6c</anchor>
      <arglist>(unsigned long int i1) const </arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Graph *</type>
      <name>graph</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a9153da43d0b9ea25b0002084f3a4568d</anchor>
      <arglist>(unsigned long int i1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a63a9de6f8d72fad03e67fe5009231840</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>erase_if_last</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>ad5e9c0d4a4c5244aa7d64507685f33f4</anchor>
      <arglist>(BP_Gen_Group *i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>erase_if_last</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a1f2695066f4133d5ab352e2094a3558e</anchor>
      <arglist>(BP_Gen_Graph *i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>aa9b8ae38b153928e364bff6003ec4a47</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Gen_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>aa0afbcdd0a98ca4b886b5fe363f8124a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Gen_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a560a7489bfb1a931f7eaa6c84b8fb2cf</anchor>
      <arglist>(BP_Gen_GVec *i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BP_Gen_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a43c80f3991c87f24a85edb24bcb9ce9e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_group_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>ae0e519f624b08a5d23cdf67d129fca21</anchor>
      <arglist>(BP_Gen_Group *addy, unsigned long int indy)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_in_group</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a89593423ef362ccf54621e91c4d0994f</anchor>
      <arglist>(BP_Gen_Group *addy)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_in_graph</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>ab049d5a98a2669df7aa83f557e50af3a</anchor>
      <arglist>(BP_Gen_Graph *addy)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Gen_GVec *</type>
      <name>BP_GVec_home</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a76ca79f13d71f34a222367c42a561e52</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>unsigned long int</type>
      <name>BP_GVec_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a02b85f2afd0a232953e14c5b6c01a7b7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_Group_Data_class &gt;</type>
      <name>groups</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a506c8fc1a1b178f5acaac7ddc69690d3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_Graph_Data_class &gt;</type>
      <name>graphs</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a516660c43c0f64b077dc1af501b6b518</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_Gen_GVec</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a8454445ca7061e03b3d0f67e4a05ff8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_GVec</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>aefad06138deb56b48f65dc4c035d8b2a</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_Gen_Group</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a01cd2ae6fc22b6d71e8c98d1539b6431</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_Group</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>ab18ab234dba7cf5a558c03b60dfb9ba7</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_Gen_Graph</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a494019c401292e716d0d8df5852b6b2e</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_Graph</name>
      <anchorfile>class_b_p_1_1_b_p___gen___g_vec___member.html</anchorfile>
      <anchor>a73997afb40a83f56b533bb95f470a131</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_GVec_Member</name>
    <filename>class_b_p_1_1_b_p___g_vec___member.html</filename>
    <templarg>T</templarg>
    <base>BP::BP_Gen_GVec_Member</base>
    <member kind="function">
      <type></type>
      <name>BP_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>ae6a6acd099278060db66fd4fefa1b32b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>abd22cd55aec4e9a89028613508ae7e24</anchor>
      <arglist>(BP_Gen_GVec *i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>a18b023141ae0b81e458ba38681bfa66f</anchor>
      <arglist>(BP_Gen_GVec *i1, unsigned long int i2, const T &amp;i3)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>abf7159ddeea2cd786eecccc31743aa6a</anchor>
      <arglist>(BP_Gen_GVec *i1, unsigned long int i2, const U &amp;i3)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>aa6c8927044c15b13b76dac38ed7c3e56</anchor>
      <arglist>(BP_Gen_GVec *i1, unsigned long int i2, T *i3)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BP_GVec_Member</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>a6ba0cb2c13ef6f71df7529e6319eb786</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>get_val</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>a55a4fc54002325e4efb9f089afc5d2f1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getp_val</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>ac807952b8788014c381141fe0b29e2f5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>alloc</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>a3b3dbaab5c0029895f7744e55c71e416</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>T *</type>
      <name>val</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>ac8625ce8952b346a6e39cceaeed31eae</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_Gen_GVec</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>a8454445ca7061e03b3d0f67e4a05ff8b</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend" protection="private">
      <type>friend class</type>
      <name>BP_GVec</name>
      <anchorfile>class_b_p_1_1_b_p___g_vec___member.html</anchorfile>
      <anchor>aefad06138deb56b48f65dc4c035d8b2a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Gen_Group</name>
    <filename>class_b_p_1_1_b_p___gen___group.html</filename>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>ada345a76d4674f56c427983296ba6d49</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned long int</type>
      <name>size_alloc</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a86be7bc99999a2ad6ca7fcddc2c26bd7</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned long int</type>
      <name>size_capacity</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a1b716c6bb658c6e1da2721a8dbe1c46c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>acf0a1615caa27709f67ea2e47dcebf7a</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>cap_increment</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a3b7c3461880410cbd11910156032c247</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual string</type>
      <name>set_Name</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a413ed97359334ba7b23f4e2db45fbb97</anchor>
      <arglist>(string i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual string</type>
      <name>get_Name</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a50e972aada29c03a7a1f781b622430bd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a7213fd1d1908d19523878878edded290</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_if_last</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a2f3a65b9c63a0881e74546b09f86691e</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>af3be4e3f1ddd94e512ed60fb04517a1c</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a60c5c7292a04236e8ca783c32f8c0aef</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a44e52cbd4ae8a3595adece061d7aa0bf</anchor>
      <arglist>(const BP_Vec&lt; BP_Gen_GVec_Member * &gt; &amp;)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>af28dbaa3f10a5d0c04b3b3e86d786cae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>affab9e02ac272e73287cf833baa75cf9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a48cf7747eeddea84aad6d555b237df35</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>ab55f23eef6893eb1e4622aa6c185696f</anchor>
      <arglist>(ostream &amp;sout)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual unsigned long int</type>
      <name>get_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a38f5c819a3e8328c53bb32c7c527dbba</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual bool</type>
      <name>contains</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>aa1ece2ff52e757881efe8bc200329aec</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Gen_Group</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a4bb1ea5f75e7d727627bd140af2cee09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_Gen_Group</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a32c6e7e2e60f6de3724610750dcf35c2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>remove_part_one</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>ad63a4de75dcff2df0e6fdead87d7110d</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove_part_two</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>ace1f8ac135ca3111fccf69fc2879733b</anchor>
      <arglist>(unsigned long int vi)</arglist>
    </member>
    <member kind="variable">
      <type>string</type>
      <name>Name</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a1df07c7098eca3fa3faef0d2b02b73c4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>BP_Vec&lt; BP_Gen_GVec_Member * &gt;</type>
      <name>group_members</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>a3ba8104b9669725d0029f526284347bb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>BP_Vec&lt; unsigned long int &gt;</type>
      <name>group_member_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___group.html</anchorfile>
      <anchor>aed99aa285d19f3fb018fd2d72ddab30a</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Group</name>
    <filename>class_b_p_1_1_b_p___group.html</filename>
    <templarg>T</templarg>
    <base>BP::BP_Gen_Group</base>
    <member kind="function" virtualness="virtual">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___group.html</anchorfile>
      <anchor>a7819beaae6d308de80d4341308d6d787</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>member</name>
      <anchorfile>class_b_p_1_1_b_p___group.html</anchorfile>
      <anchor>a490454dc83025789e5e5505163fe2de4</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual BP_Vec&lt; BP_GVec_Member&lt; T &gt; * &gt;</type>
      <name>all_members</name>
      <anchorfile>class_b_p_1_1_b_p___group.html</anchorfile>
      <anchor>a3c04b2b854f60ae0efab8beeb302f784</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual T *</type>
      <name>getp</name>
      <anchorfile>class_b_p_1_1_b_p___group.html</anchorfile>
      <anchor>a0775e7aa190f909ba1d64c277f14c1a1</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___group.html</anchorfile>
      <anchor>ad10b74b30779f73f2507622ee8ab81c6</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual const T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___group.html</anchorfile>
      <anchor>aac1ed0a0ea5a6b07cdc64b049e5624db</anchor>
      <arglist>(unsigned long int i1) const </arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Edge_Data_class</name>
    <filename>class_b_p_1_1_b_p___edge___data__class.html</filename>
    <member kind="variable">
      <type>BP_Gen_Edge *</type>
      <name>address</name>
      <anchorfile>class_b_p_1_1_b_p___edge___data__class.html</anchorfile>
      <anchor>ad19d0fec40fc91aeeb6793aaa07ed9dd</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>unsigned long int</type>
      <name>vert_index</name>
      <anchorfile>class_b_p_1_1_b_p___edge___data__class.html</anchorfile>
      <anchor>acd4fcd95988c238e5fe2f580f08c1c6d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Gen_Vertex</name>
    <filename>class_b_p_1_1_b_p___gen___vertex.html</filename>
    <member kind="variable">
      <type>unsigned long int</type>
      <name>vertex_list_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___vertex.html</anchorfile>
      <anchor>aa52c74216ede835445535d126878c537</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Gen_GVec_Member *</type>
      <name>val</name>
      <anchorfile>class_b_p_1_1_b_p___gen___vertex.html</anchorfile>
      <anchor>adf15dc7f16d1aa0e261b8a10c5eba7d9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Vec&lt; BP_Edge_Data_class &gt;</type>
      <name>edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___vertex.html</anchorfile>
      <anchor>af84bbcf7285bdd6d2cd51516b30018b4</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Gen_Edge</name>
    <filename>class_b_p_1_1_b_p___gen___edge.html</filename>
    <member kind="variable">
      <type>unsigned long int</type>
      <name>edge_list_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___edge.html</anchorfile>
      <anchor>a5cbb41de3d12220a0d1c243b03e84430</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>bool</type>
      <name>dir</name>
      <anchorfile>class_b_p_1_1_b_p___gen___edge.html</anchorfile>
      <anchor>abc712f05b4e69b345e112879c26f93c7</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Gen_Vertex *</type>
      <name>vert_list</name>
      <anchorfile>class_b_p_1_1_b_p___gen___edge.html</anchorfile>
      <anchor>a190d2cd97ea8b86a041072f361aec85d</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable">
      <type>unsigned long int</type>
      <name>edge_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___edge.html</anchorfile>
      <anchor>a4a3a8d4594cd1d8dbbcfae291afbba61</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable">
      <type>BP_Gen_GVec_Member *</type>
      <name>val</name>
      <anchorfile>class_b_p_1_1_b_p___gen___edge.html</anchorfile>
      <anchor>aa48c971a374609fb4068b08e43108034</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Gen_Graph</name>
    <filename>class_b_p_1_1_b_p___gen___graph.html</filename>
    <member kind="function">
      <type>unsigned long int</type>
      <name>vert_size</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ac6ed317aa17a771d8b53c47d081c9965</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>edge_size</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a960ad1aaffdae4c2ad16718ca33239f6</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>num_incident_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a2ed38775eb015dab685b1d70c1d89671</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>num_incident_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a45d4fbf8f4470bc1b73383f061a74709</anchor>
      <arglist>(unsigned long int i1, int i2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>num_incident_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a93241a87263c91be72539522f2b0d85d</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>num_incident_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a6792bf839ece36e5881fe673e7f17ab3</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, int i2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>num_incident_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ad3c92a31823e9a887dbed55ff7bd484f</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>num_incident_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ab1247adb0576bcd9fe3664ec58779059</anchor>
      <arglist>(BP_Gen_Vertex *i1, int i2)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual BP_Gen_Vertex *</type>
      <name>add_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a455dbc6ebcb4f72eec10da6cbdfd2d36</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a262bbe5e03d2de26db5c61c0d0e96f8e</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_vert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a72ad6644d3b10b16840776cc19d484d1</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a2c10461e479225320201d6b7194e3c58</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a814d87e25dd547072cb7163eab482ad8</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>acebe5911793711efb4cc36608077213f</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a4bb9d45ce3d2c1ee9bff1b32d2e16052</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_if_last_vert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a639280cfe3310ecaa93435188827f955</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_if_last_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>af110af14168f43792fbdded13f5e8f95</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_if_last</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a2f9fb6150cad1632394700f8f50bfd61</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_if_last</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a80d7aeb2d6dbdfa51eb26647a4a5879e</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>erase_if_last</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a7f7e35f1a682a10ad5bfb82ca8f3fe19</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>aa7ea888bb04b411145d56d8b94dc1eb8</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a5394784a2839d78f39ab9664bf872bc6</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ae31af1e85452573bae86c6e1e92bf734</anchor>
      <arglist>(BP_Gen_Vertex *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>disconnect_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a9e294687ed2aca0444fab891871db968</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>disconnect_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>aed2eca5b53110eaea8b28deb6fb01c1d</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>disconnect_vertex</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ac3ac1a66f24e4925fa4b92b1c6ef1760</anchor>
      <arglist>(BP_Gen_Vertex *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>add_halfedge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a86fcb9ac14416327dc2c8847cb1d3da6</anchor>
      <arglist>(unsigned long int, BP_Gen_GVec_Member *, int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>add_halfedge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a4e4d99bfd86268610a35042e8f6275ca</anchor>
      <arglist>(BP_Gen_GVec_Member *, BP_Gen_GVec_Member *, int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>add_halfedge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a4274d3f995e6848cc29f29c298fb8fc7</anchor>
      <arglist>(BP_Gen_Vertex *, BP_Gen_GVec_Member *, int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>connect_vertices</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ab69c6eca5767200c4301800eca7d5aec</anchor>
      <arglist>(unsigned long int, unsigned long int, BP_Gen_GVec_Member *, int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>connect_vertices</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a7d1dc6a85b99fd3288740c051f06fc2b</anchor>
      <arglist>(BP_Gen_GVec_Member *, BP_Gen_GVec_Member *, BP_Gen_GVec_Member *, int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>connect_vertices</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>afbd9158566bf7649671020a9ea514cac</anchor>
      <arglist>(BP_Gen_Vertex *, BP_Gen_Vertex *, BP_Gen_GVec_Member *, int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>disconnect_vertices</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a8e4e0827f52090e058e4f54793578f41</anchor>
      <arglist>(unsigned long int, unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>disconnect_vertices</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a107734137a6d36ac40b2116625a87dec</anchor>
      <arglist>(BP_Gen_GVec_Member *, BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>disconnect_vertices</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a5f36ce4b8c7e685700c81d0478d283fa</anchor>
      <arglist>(BP_Gen_Vertex *, BP_Gen_Vertex *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a667a3206ce7485bce790e2c4e5fa7597</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ad1928e52881ecf7e8a6244b05ce72dc0</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>remove_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ac8f811868dc2389aa98259c4973b1c64</anchor>
      <arglist>(BP_Gen_Edge *)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clear_edges</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>aeda420d6f7bfffde6cace22ec5cb8623</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ada8f75a455cf27e770dcbd8b031a33b2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a22eb7278bc1842dd28ff43b8748ade8c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ac6aad20d45e058c5357555b5567aee4e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a47aa198022e578bef40ac10fc4cdfd36</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>cap_increment</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a5375011970f1d1e302bbdf0518ccfa01</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>vert_to_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a83e5b6eb2386bad0ee20c066b5f540cc</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>edge_to_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a43cd1a3b60813275aabe83515033acf5</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Vertex *</type>
      <name>index_to_vert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a9494d83fd274753b906a069845f9eb96</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Edge *</type>
      <name>index_to_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a0d19c53569a3e2f80e5e3d894046cb78</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>member_vert_to_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ab6a504ad03a01d9dfa2d38df6ec22c3b</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>member_edge_to_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a5f0dec714a9c489d5d3bb67daef3ec4a</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Vertex *</type>
      <name>member_vert_to_vert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a8f024c977e95d1252825ee352139419c</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Edge *</type>
      <name>member_edge_to_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a436d91b6f6fd6c1180d2b6ebb55dcc8b</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_neighbor</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a9104ff121afdd5cec565f6546935a264</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_neighbor</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>af83385896e17b787efd0275faaa893ea</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2, int i3)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Vertex *</type>
      <name>get_neighbor</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ada6bc77e4cc4731ee7741adbc541370e</anchor>
      <arglist>(BP_Gen_Vertex *i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Vertex *</type>
      <name>get_neighbor</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a0b62b0317cf1e80574148a6828aee357</anchor>
      <arglist>(BP_Gen_Vertex *i1, unsigned long int i2, int i3)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_incident_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a01eb142d4509665edec3e9d0d83cd060</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_incident_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>afb7feaa1694d13ea82e6b767b0883b4d</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2, int i3)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Edge *</type>
      <name>get_incident_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a009a76f23c676d3102fad6a7b42e2d52</anchor>
      <arglist>(BP_Gen_Vertex *i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Edge *</type>
      <name>get_incident_edge</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ae366ae6ff846350bd75f848d46b59a44</anchor>
      <arglist>(BP_Gen_Vertex *i1, unsigned long int i2, int i3)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>are_1NN</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a999152e5b56393d8def4189c8389ade5</anchor>
      <arglist>(BP_Gen_Vertex *i1, BP_Gen_Vertex *i2)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>are_1NN</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a22d33ea55895bdeec7c384b66afab61a</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>are_1NN</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ac2152b4298826dded4b78ab60453cd60</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, BP_Gen_GVec_Member *i2)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>edge_startvert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a9e14335c920662205fcbc712417b295c</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Vertex *</type>
      <name>edge_startvert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>ab33f96d943f05bc47a4a03ab2cf6155c</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>edge_endvert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>aee58cfb885ea14398b27b6d31445c3be</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_Gen_Vertex *</type>
      <name>edge_endvert</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a8d0c58f865ca1899462abf5d5541530a</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>contains</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a0b01d10e4af1cea2787d8681542669ff</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>aea1b64d81814fa928a3b8ba26b597b3b</anchor>
      <arglist>(ostream &amp;sout)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Gen_Graph</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>aa1a818e038953964de62e6a259902f08</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_Gen_Graph</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a87cdfdc16d6599a12f46c16486b1a346</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable">
      <type>string</type>
      <name>Name</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>a29e73f8c6676facbc87fbdc6f5b6363a</anchor>
      <arglist></arglist>
    </member>
    <member kind="function" protection="protected">
      <type>unsigned long int</type>
      <name>get_index</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>af0423dc93febdf417ec048b0b3bbf997</anchor>
      <arglist>(BP_Gen_GVec_Member *)</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>BP_Vec&lt; BP_Gen_Vertex &gt;</type>
      <name>Vertex_list</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>abe14a941873b0d0f16b7eae5359a65cc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>BP_Vec&lt; BP_Gen_Edge &gt;</type>
      <name>Edge_list</name>
      <anchorfile>class_b_p_1_1_b_p___gen___graph.html</anchorfile>
      <anchor>af469fad7d7c453b1110b069009d50073</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Graph</name>
    <filename>class_b_p_1_1_b_p___graph.html</filename>
    <templarg>T</templarg>
    <templarg>U</templarg>
    <base>BP::BP_Gen_Graph</base>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>index_to_member_vert</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a61276afcfc9e6071ff5924e4ab9e30b9</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; U &gt; *</type>
      <name>index_to_member_edge</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>ac3a98c88981e45c7b4eddcd8c84c618c</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>vert_to_member_vert</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>ad2e183e9b65c037e960174ec9239e58d</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; U &gt; *</type>
      <name>edge_to_member_edge</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a1a8a06910717c6a036145995e73ade81</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>edge_startvert</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a878986a79a3dd23235cd92dea299699d</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>edge_endvert</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>ac653e6ac17829a02e3e49f2da9d7d5b1</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>vert_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a6e306e13415bc7a468ccf804f2e7c22f</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>U &amp;</type>
      <name>edge_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a30a413752cea278cb24eea569335477d</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>vert_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a633c5f2f92b07bc4b7ff738d43d1aac7</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function">
      <type>U &amp;</type>
      <name>edge_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a1a6f93d664fa390c4b87e78c4dddd36a</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>vert_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a4de3a090a2af543bff44ed24835d6edf</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>U &amp;</type>
      <name>edge_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a9e234b82fcc9949e9a8b6fbbb2fcc0af</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getp_vert_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a4aeafcd8b448254a3229cbbea007cdb2</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>U *</type>
      <name>getp_edge_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a7039924f0e2c930af8404ae50e6db558</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getp_vert_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a5f78ced445e1cf82cb842cfefae8ada3</anchor>
      <arglist>(BP_Gen_Vertex *i1)</arglist>
    </member>
    <member kind="function">
      <type>U *</type>
      <name>getp_edge_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>abddaeb5881d1f7562b4bd0a7382c9ffa</anchor>
      <arglist>(BP_Gen_Edge *i1)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>getp_vert_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a75036d6b831739a6d9bd172e6c20c7b4</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>U *</type>
      <name>getp_edge_val</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>aaf259c26aba4f8be4f781e6decceb743</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>get_neighbor</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>ab6a620979e1360cb52da82e3c8274f3d</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>get_neighbor</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a288156072a75ca693b1f944f45542221</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, unsigned long int i2, int i3)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; U &gt; *</type>
      <name>get_incident_edge</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a249b373b6799921f1dc7a6e9b72d14fb</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; U &gt; *</type>
      <name>get_incident_edge</name>
      <anchorfile>class_b_p_1_1_b_p___graph.html</anchorfile>
      <anchor>a6c846c7b83281891c6d8cf384d16df70</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, unsigned long int i2, int i3)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Parse</name>
    <filename>class_b_p_1_1_b_p___parse.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_Parse</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>ad7364df139b2039966f1f2611cb46b7e</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_Parse</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a877f8ce3d412567cca8d116a5121eb6b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a418e47eb4b2a802785613fa773028bc1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>af08c7a1d9265d0e3843eda599291495a</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_com</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>ae42b78c4a8796385c3a45290c1ff7e83</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear_com</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a7be00c29740f297040ec175194138974</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove_com</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a8bc4ffa172aedf3b9a008ebedae6ac7c</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>get_com</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a553fb9b1d6f3481410185ee98e9fe79f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>com</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>aa91dc9fed05592a45f456a8b9627ed4b</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getline</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a9e51134d6d6dedd0b752ea6b1490ca47</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getline_all</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a789d6f79a12b9783cc65f8101bc9cacd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>getline_string</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a65b11df7740886eac25cfcf89a9ff488</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>getline_int</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>aefbe2bc0b5849dbffd35c7c79e056722</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; double &gt;</type>
      <name>getline_double</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a4a10785352dccfeeeff02c9d22adb6e9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>next_string</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a8675a87991cca4f083934336df4415c0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>next_int</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a5fcfea06fa4ba2f48fb2727e481536b4</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>long int</type>
      <name>next_lint</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a744c4f1f7e26cf0599f48fff63a29b54</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>next_ulint</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>af2e1d6c56fbb8c2ae03d5704a10e6bf2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>next_double</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a639a655a8fa0e64239be38baf3aba1ba</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>eof</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>adad55c254b4dd5e8b8eaf0d48f1b772f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ifstream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>a03368759b37e1951586a4fcd6142b6eb</anchor>
      <arglist>(T &amp;t)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ifstream</type>
      <name>infile</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>abd16ad2c21b4e7b676d6761aaf2ff8e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; string &gt;</type>
      <name>com_list</name>
      <anchorfile>class_b_p_1_1_b_p___parse.html</anchorfile>
      <anchor>acb6f9c7a44b94e0c5b9c393ce5163014</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_bParse</name>
    <filename>class_b_p_1_1_b_p__b_parse.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_bParse</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a2f2ca00c7e82388b412647976c3570c7</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_bParse</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a9c7aef1eb18140538831d643cf9ca093</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a9896dde8914dec02c09d027122ad8939</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a19c95c26d0f298642f00d4ce8410c7a8</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type>char</type>
      <name>next_char</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a8cb4f823d3ada302406ae7865e985b3d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>next_int</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a6d85713386b7f68d0c48246feac6fc5e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>next_ulint</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a8f1e29ff3ef74d8d3171b91d81c0deeb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>long int</type>
      <name>next_lint</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a7cf948e204759d737221cf700f8de17d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>next_double</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a21e6bf63f8f0187ce3dc0432ca8f76cd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>peek</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a0025da152fbb4f78a4457f4eedd6a5e5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>eof</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a1c28734960eb1d083b242a3f67e0849b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ifstream</type>
      <name>infile</name>
      <anchorfile>class_b_p_1_1_b_p__b_parse.html</anchorfile>
      <anchor>a30797cfeba657fffe4c1356c3658d88e</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Write</name>
    <filename>class_b_p_1_1_b_p___write.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_Write</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>af136ad1690f0c72e2fb63fe3fc66737e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Write</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>ace0d36cffa4da6cc50d616a6e7d2cdb3</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_Write</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>af6b235a96e50e7a3e8ec08eab2ac736d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>af186a44682cf7705f8a0f9c4bfcf4209</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>a5e2a6473412c0680b68f0cae230d8233</anchor>
      <arglist>(string s, string mode)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>newfile</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>ad8f6899389ef96a686a795a569035b8b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>name</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>ab9bd589116f7c1bda03e1f1f291bbad0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>ad2a82aae09b39e580e1ae405a8ca00f9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ofstream &amp;</type>
      <name>get_ostream</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>a9ba6ecaa0ba267775e663f11fccb2e50</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ofstream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>a86c44bed2076a9844331f73deac1d995</anchor>
      <arglist>(const T &amp;t)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>append</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>a42857f2cf1fabe05cc3a3a092302633f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ofstream</type>
      <name>outfile</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>a569440716608bf5fd8bd47b56c9f22e9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>filename</name>
      <anchorfile>class_b_p_1_1_b_p___write.html</anchorfile>
      <anchor>a7423b22fcb99648c79ae13ff904a1e39</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_bWrite</name>
    <filename>class_b_p_1_1_b_p__b_write.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_bWrite</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>aac756bd28835b819d0dfa403f42a455f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_bWrite</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a2eda968dfad25f2a10872597b4a5a70e</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_bWrite</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>af5fce244846536975d146d5d7f95cadb</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>ae69569d89dcc3e51219397ae427b8b0d</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a5036d5a676d28f38fbd8a925422ca952</anchor>
      <arglist>(string s, string mode)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>newfile</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>ab503da79b669b46514e3040a69d0cdd2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>name</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a216b7492f33ef9841cd69a762a73cf84</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a7cd021fc39e20301c7b1e00e001b4311</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>ae435c0a2c5adc508d17210c4b6326c1d</anchor>
      <arglist>(const char *s, unsigned long int n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_char</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>ab6db15f8249fa114b9ddb43ae991423f</anchor>
      <arglist>(char c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_int</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a6f9a39c3ab4db158e024f6608db4300b</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_lint</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a13c5122a9bafe4866600e8edfffe501f</anchor>
      <arglist>(long int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_ulint</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>aa87c41f62e6b543ebb348f8a1af2eb9f</anchor>
      <arglist>(unsigned long int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_double</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>acf976dea472ac968328c14b25bdd7ff7</anchor>
      <arglist>(double d)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>append</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>ab17e84dbb621e6b4cc02c3a568ce1f0a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ofstream</type>
      <name>outfile</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a1a7279354dfdb2fde65f391807271413</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>filename</name>
      <anchorfile>class_b_p_1_1_b_p__b_write.html</anchorfile>
      <anchor>a2cc4eea51bde9d418efeb9a0450c0dec</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_RGB</name>
    <filename>class_b_p_1_1_b_p___r_g_b.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_RGB</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a705282104843585946c15af6c59d5541</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_RGB</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>ab13ea568cf2a523944f7d7d1d323f4d8</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_RGB</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>ab7faa60296a6f9a3ae8a45d1f578e399</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_RGB</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>aeb950d9385805912880b2ec73280c114</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a1f5d2eedad608df659761ff1f5c0a36e</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a234e8148f7221258a944372a25ab9bda</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a8f16a38635b08d3dc6920a633d4a3d38</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>double &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a65ec561ca48f60121c5e48f44bb97dcf</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>const double &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>af641bb5d92824eb5604d99a1f35ca3b3</anchor>
      <arglist>(int i1) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_color_from_map</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a9d1d001b7d79629e1a5145e5274ebc5c</anchor>
      <arglist>(double d, string map)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>val</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>a7f964e460332dc261614ade027124377</anchor>
      <arglist>[3]</arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1_b_p___r_g_b.html</anchorfile>
      <anchor>ad7cea6b28dc639ba6a3928f553d0ea67</anchor>
      <arglist>(ostream &amp;outstream, const BP_RGB &amp;col)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Text</name>
    <filename>class_b_p_1_1_b_p___text.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_Text</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>aeb2803f12f04ed4727cfd70ae54d0d92</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Text</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>ad91253707649858a5c176d14715d1b27</anchor>
      <arglist>(string s_text, string s_font, double d_fontsize)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a48cee2e6c7248dc9dfeb78ff97f60da1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>af09942ef1d080818b5fe1b827be678bc</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_format</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a9697376a2fe79ddd59afa4f23a5a7a47</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_text</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a0418dd0ca53886312047cd342af1f79d</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pos</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>abd15f1a57e4124e6acad6f808218bf74</anchor>
      <arglist>(double d_pos0, double d_pos1, int i_pos_mode, int i_align, double d_angle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pos</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a8327736cfb91e0bad42880bf7d3ba08b</anchor>
      <arglist>(double d_pos0, double d_pos1, int i_pos_mode, string s_align, double d_angle)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fontname</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a13f38e66efe8de7a116ca25336847a6c</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fontsize</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a015b692fd4029284db2687d91df51493</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_text</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>aec9d88753fa70ca7ef30eaf3a9d2e2e4</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_color</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a5ae4decfc776753aeafa34de5ed1808d</anchor>
      <arglist>(BP_RGB c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_color</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>abe3055b27feed179dc3beb60b8001a2d</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_color</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>acf1c7225010e0baafe014bc77c6f0572</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_color</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a67a6c4528a2be929082b96e1682e8047</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_textpos_mode</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a6b78bdef8a1858725aac61a3ce04aa77</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_textpos</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>aeb597dd13c1c1424de45e955f9e586ad</anchor>
      <arglist>(double d1, double d2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_text_angle</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a469dfb7c3e6afc77f735e1af32558cf3</anchor>
      <arglist>(double d1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_text_align</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>ad0492dbe5a2948db11e5d53ccffcc80d</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_text_align</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a5d85baf65c064897f01d07ce3c2496e3</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>text</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>af8ecf43273abe3487e99ac45f6e3f9dc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>fontname</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a9848032ab0ba18c4f3ae37e204d04bc0</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>fontsize</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a6b0d3f87f191584ce11f6f0a96b38210</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_RGB</type>
      <name>textcolor</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>ada1b0449247ff64a35ec1a5d5d2552fc</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>textpos</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a9925b151ba4cf306bad8154c93e0db24</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>textpos_mode</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>aa53393116cf09f6e1e2c3c18557f6f5a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>text_align</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>a72072745a551e878e1c2e5c6e05b19a4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>text_angle</name>
      <anchorfile>class_b_p_1_1_b_p___text.html</anchorfile>
      <anchor>ad77d9cb35a9c46b40672280b4df64044</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Plot_Data</name>
    <filename>class_b_p_1_1_b_p___plot___data.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_Plot_Data</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ab25b606d96df981c6f30d86d772ffdbf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a46c22ed8a5e7091effe65a4612f5f544</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_ps_data</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a078108dc582539fbe7abb4bb6caea7df</anchor>
      <arglist>(double *axis_pos, double *axis_size, double *axis_lim)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>aa6e64cbdd12167e8693d459759f39fe0</anchor>
      <arglist>(BP_Write &amp;file, double *axis_pos, double *axis_size, double *axis_lim)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_data</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ae43c211195725f158ee7656668cd9de2</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdat, const BP_Vec&lt; double &gt; &amp;ydat, const string &amp;nam)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lineprops</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>afdd2a08b2ba63f10aad180164b449341</anchor>
      <arglist>(int color, int style, double width)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lineprops</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>afdc3cc83aacd9a3700cfce8ef60324fb</anchor>
      <arglist>(string color, int style, double width)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lineprops</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a42da3e645ea8a84969c2b11508049bd6</anchor>
      <arglist>(BP_RGB color, int style, double width)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_line</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a606dfe42fb8c302f186c6a270da49fc5</anchor>
      <arglist>(bool i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_linestyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a5c3cd76c7754d418f6089d6292b17f79</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_custom_linestyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a457f44b8abedfc3c653dccd8937631f4</anchor>
      <arglist>(BP_Vec&lt; double &gt; &amp;i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a73bed5dc362de53751a421689bad5064</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_linecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ab86d291f51fa2a6f6b6fe5f4c2f2921e</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_linecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a68f0ceec421e36d21f342f948c7099a9</anchor>
      <arglist>(string i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_linecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>aa1951e17d1c810f99c4ad90da3f8c096</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_linecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ab58dd09cc1d3b2d1d2e4f71e86881059</anchor>
      <arglist>(BP_RGB i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointprops</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a8bdbcb0ee98366c2306d07046fb661c9</anchor>
      <arglist>(int color, int style, double psize, double pewidth)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointprops</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a8dec25410809e27eff66446062336600</anchor>
      <arglist>(string color, int style, double psize, double pewidth)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointprops</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a6c1234805a1d6c7602de775bbd7d133d</anchor>
      <arglist>(BP_RGB color, int style, double psize, double pewidth)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_points</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ac2d9042841b9533283924f129a7d1e79</anchor>
      <arglist>(bool i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointface</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>afb85b290ee98460f1d76a88e7264961e</anchor>
      <arglist>(bool i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedge</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>aee169fcb2eda0b60ed3a44d126c81c2e</anchor>
      <arglist>(bool i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointstyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a2d04282a19f8dcbffa9a9697306a80cd</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointsize</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ac3574f09ab277150fb3cf6a9b5b5dbd5</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgewidth</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a11c4983e5a8fa71e4e1963c82bc7aea8</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointfacecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>aa3f94d9af5d6e36e761ddce5e85efb64</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointfacecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a780c018478762b7d8a88ecef9e06cfe7</anchor>
      <arglist>(string i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointfacecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>abaa0d7629809db3fafca2376a0a1abe8</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointfacecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a1501b3f5b1be5123215960dd817b2d6d</anchor>
      <arglist>(BP_RGB i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a92f47629ceeff82be25f3b6d30f1b619</anchor>
      <arglist>(double i0, double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a2f92a28543be3d52c07e1b642f0335d1</anchor>
      <arglist>(string i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>afa97b8b743224031078455cd8d058ff0</anchor>
      <arglist>(int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a7e9e938d189dde3ef8026e6983569c29</anchor>
      <arglist>(BP_RGB i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointstyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a9c94ac00ee4c5cd5c6fc14599cb339cf</anchor>
      <arglist>(BP_Vec&lt; int &gt; i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointsize</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>aadc211adfe08388704c101592e79bb0a</anchor>
      <arglist>(BP_Vec&lt; double &gt; i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgewidth</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a6f9a9b88b3983e0d05ec3a65f8946f2b</anchor>
      <arglist>(BP_Vec&lt; double &gt; i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointfacecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a959b6b6f9739b9f6b08e1652263ac6fb</anchor>
      <arglist>(BP_Vec&lt; BP_RGB &gt; i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_pointedgecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a802e5c7a4e20a1f0f37578ea2d3243ac</anchor>
      <arglist>(BP_Vec&lt; BP_RGB &gt; i1)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>name</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ad74c141e129358f9920824ff3a727f1a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>xdata</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a4a5b8d543c5d7e64324a528333b40a4e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ydata</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>afe69e3a4447951cb269a74539ed0bcfb</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ps_xdata</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a4810e41e42dc01004d269c51a7b8dee9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ps_ydata</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ad671ca816c93f9fc00344768f86c948e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>line</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a89a3c3132b5c87e6f7110a137d214ee2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>linestyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a852c8073646f27f3cdea5b1a5a24e82c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>custom_linestyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a7b3b93059e77f35c9ce4c45d520b8ee9</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>afb463fa803d678a96e2b283d44fb4d2c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_RGB</type>
      <name>linecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a7a271a80945879936eb6661dd5b0d8c5</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>points</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a44545721038348c8f2ba859420293333</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; bool &gt;</type>
      <name>pointface</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a5230db156e9435afcf98d95aaf4bcd60</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; bool &gt;</type>
      <name>pointedge</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a1381d55a3f1430a5ab81de032c3a975a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; int &gt;</type>
      <name>pointstyle</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>ad2b4d2710dcab0747981d15acd122494</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>pointsize</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a69e9599992aec35b8ea148def1eec914</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_RGB &gt;</type>
      <name>pointfacecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>a518dc9bd6c3f15fe2d7ec69709397de1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_RGB &gt;</type>
      <name>pointedgecolor</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>aa9c14e60fbb07d05c7da120a366f8707</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>pointedgewidth</name>
      <anchorfile>class_b_p_1_1_b_p___plot___data.html</anchorfile>
      <anchor>af192ed436cf3fea719c0eeec362b1acb</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Axes</name>
    <filename>class_b_p_1_1_b_p___axes.html</filename>
    <base>BP_Vec&lt; BP_Plot_Data &gt;</base>
    <member kind="function">
      <type>void</type>
      <name>set_tick_ps</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a2c36027dfb29a1a9f31ec806913eaf09</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Axes</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a40c9a7d0069d3b1e894ca9b02734c69c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fontname</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aeb59131a565bd806007beacdffad0e8c</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_fontsize</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aaefddce18b995da8c5b753198e239e50</anchor>
      <arglist>(double d)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ac232f0640293bfa0371c85871a58220c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_auto_ticks</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a8cf1f6ca6f5bcf68844e868a4b1c2191</anchor>
      <arglist>(int axis, double lo, double hi, BP_Vec&lt; double &gt; &amp;tick_val_major, BP_Vec&lt; double &gt; &amp;tick_val_minor)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lims_and_ticks</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a76f754eddb4c6617a54a5738954dd371</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_background</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a2b290412d9edafdba0f9d8c04ef766a1</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_edges</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a64a8ae4152d391df071fe5d26f6b6e65</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_ticks_and_ticklabels</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a875a12760e150daebceb3e0a94468447</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_axis_labels</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a0d44fd7361cb3c214d64e5c6fb74453e</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_data</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a65c7de877bdc26af10a00df6419b71ea</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_legend</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a1d6a2bf1cd724956bec0d0775b774820</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ada9ed9dae8ea2783a791013524899137</anchor>
      <arglist>(BP_Write &amp;file)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_data</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ab62c3c15817e090a29ad4566d1c9bd19</anchor>
      <arglist>(const BP_Plot_Data &amp;data)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_line_wpoints</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a0c3a390ea4840d5f327a61299bd418f5</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdata, const BP_Vec&lt; double &gt; &amp;ydata)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_line_wpoints</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ae0ac5398f619b4d2818f7fb10cae4439</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdata, const BP_Vec&lt; double &gt; &amp;ydata, string name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_line</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a26778d35d79c68f912e4c29127c88e75</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdata, const BP_Vec&lt; double &gt; &amp;ydata)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_line</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a0f2a683083b0a62a74fef9b8fbb34542</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdata, const BP_Vec&lt; double &gt; &amp;ydata, string name)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_points</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a77e00c1dd669a403684c12250a94cf9e</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdata, const BP_Vec&lt; double &gt; &amp;ydata)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_points</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>add55399230b3c9f0e96783cfa001302c</anchor>
      <arglist>(const BP_Vec&lt; double &gt; &amp;xdata, const BP_Vec&lt; double &gt; &amp;ydata, string name)</arglist>
    </member>
    <member kind="function">
      <type>BP_Plot_Data &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a923f79af64bef8ee2444bb8d19291477</anchor>
      <arglist>(const string &amp;name)</arglist>
    </member>
    <member kind="function">
      <type>const BP_Plot_Data &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a0160c8466b24c8d49612f72fefdef443</anchor>
      <arglist>(const string &amp;name) const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_auto_tick</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a0e067928cd5f9e6ff3e5a0a1bf13b5a3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_auto_xtick</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ade0f3d7cc55864b7bd2737bba3c2eb75</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_auto_ytick</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ac4c2de4f4db2b655318cb484c3a088e9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_xtick_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>af556bb424153a9d82b6b605df248fb43</anchor>
      <arglist>(BP_Vec&lt; double &gt; i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_ytick_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aac940a9458025e3eb9a48c113c09cd84</anchor>
      <arglist>(BP_Vec&lt; double &gt; i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_auto_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a9a19f60890dffafdb78cacc9f58d8ba0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_auto_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a4d648153121926e3589f001760c04e8d</anchor>
      <arglist>(bool i1, bool i2, bool i3, bool i4)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ae1a6b6d423d6bc4f45a434f582e1af4e</anchor>
      <arglist>(int i, double d)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a8d3fa7ccc3751e5a49b7b05e302a9349</anchor>
      <arglist>(double i0, double i1, double i2, double i3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_xaxis_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>afd66eec8e7145dd3cd5a62d61745b12b</anchor>
      <arglist>(double i0, double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_yaxis_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a09f76e1b15c7effa5c89ca6689c3ad01</anchor>
      <arglist>(double i2, double i3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_data_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a52dd465f42400edd8c936f0d748fa257</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_data_pointsize</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a744a7a982fc8a8481df11009f68f020d</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_data_pointedgewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a81c528aef10e5d746d03539c51bba55a</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_size</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a6b044a572927f03744fb6aed98fb4f14</anchor>
      <arglist>(double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_pos</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a8dc32006e07d4b6f9c662f77ba5a341c</anchor>
      <arglist>(double i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_axis_size</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a240c173a97b22ddbcc935ceb54bc7aa9</anchor>
      <arglist>(double &amp;i1, double &amp;i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_axis_pos</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ac7c14ef6a6f9239775f789e798c53a10</anchor>
      <arglist>(double &amp;i1, double &amp;i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_edge_on</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aba5e68f0a82759b2a70c617d97fb5cae</anchor>
      <arglist>(bool i0, bool i1, bool i2, bool i3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_tick_on</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a1d4f0ff1ad0866a486c51d7478044901</anchor>
      <arglist>(bool i0, bool i1, bool i2, bool i3)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_target_tick_major_sep</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a58597b0b4cfbd966e67e13c29c229a66</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_xtick_pos_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a79350d1406da6b68f00eb0d19652d774</anchor>
      <arglist>(double i0, double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_ytick_pos_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a68434d4a741c9488f9fda5a9a7237049</anchor>
      <arglist>(double i0, double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_xtick_pos_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a9cf9b516cca75273417eb992f32673be</anchor>
      <arglist>(double i0, double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_ytick_pos_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a6ff0f2bb6e9f071feb5750141ff356fe</anchor>
      <arglist>(double i0, double i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_axis_edge_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aa0d8117de3f81a508aab5cd4bac46749</anchor>
      <arglist>(double i0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_tick_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aeeb7a2995cd0eabdc83be7ab2d01922f</anchor>
      <arglist>(double i0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_tick_major_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ab5b4f461da3d7b6e7b58c0761fb29de8</anchor>
      <arglist>(double i0)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_tick_minor_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aa81af3579c5e1757c716eb001c7e7505</anchor>
      <arglist>(double i0)</arglist>
    </member>
    <member kind="variable">
      <type>BP_Text</type>
      <name>xticklabel</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a562e8ad14b284e13e016f9348861df0d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Text</type>
      <name>yticklabel</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aa72fe1e5f148c8b157ca53355adb746d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Text</type>
      <name>xaxislabel</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a44fcbaeb0b2946249687b1c60e68f6ab</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Text</type>
      <name>yaxislabel</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aff05b67090251ac6525f2cc39d14cc0b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable">
      <type>BP_Text</type>
      <name>title</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a755415e4e589cc134ec2e374bd358457</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>axis_size</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a6499413dcee4e4b9249cd08b74325e1c</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>axis_pos</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a81a5165df6906749ccc21fb3d8d9b76f</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>auto_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>af48255ce8c97ed9aedb9296b64d21a9b</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>axis_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a747e3c0d9bd59aac0409bd9fb8fcc05a</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>data_lim</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aeef8f3fa3ba05fde5bd5edc704a1fe05</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>axis_edge_on</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ac594cda044f912daa2ed3413d606b48e</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>axis_tick_on</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aa18c5f7935615c0b62447a4bd12a3a2c</anchor>
      <arglist>[4]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>auto_tick</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a3a54c113e51b64bdae41890c7c5727f9</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>xtick_val_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a00290af61afe06e574f93604f3d0e24f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ytick_val_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a3e963579097f2b0c80f1d96d2bc23fb6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>xtick_val_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a62724743ee0becbb2aad4b168323d14b</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ytick_val_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a595bdd02b8266a0bc5602e4b75af0ae4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>xtick_pos_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a5f1b88f9bb3ce482335a09c34b2afbfa</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>ytick_pos_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a8c55e350a442dbe2ff042aefc3be36e4</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>xtick_pos_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ad4f612f52f1f86200a7965fffc3c8d75</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>ytick_pos_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>acf77261ef03ad965ad40a63bbe8b09a0</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>axis_edge_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a3f8c6e35ff97c842fc208188676be54e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>xtick_major_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a64eebc5c6820d1d31c5bc5378e6ba887</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>ytick_major_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a15e5ee0ecc66940cf91561b44584131a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>xtick_minor_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a9c92d3649ca85c5037889bc9b05eefec</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>ytick_minor_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a41d17af36dcc4bfc4c470ad1ae52988e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>target_tick_major_sep</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>aeadf03def668b2c62caf843d65db5eb4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>data_linewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a66d65b36e58644fd89c51afd2cadc4ec</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>data_pointsize</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a6862e474cb7a731b283236dd6d9f221a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>data_pointedgewidth</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a300d65ab99bef3bb0ba4ea3317f525c6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>xtick_ps_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ae06dc510a21e58d1a8905429f9e8d123</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ytick_ps_major</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ad21cd2eb293161b888982f3a9c0c46b6</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>xtick_ps_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>adbbec6b0ce0d265df87ec4f73807a39c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>ytick_ps_minor</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a00ad38ba08f75cc99afcbe331bdaba01</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>fontname</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>a3a5ac70e15ff30976134f1dfeb8dc193</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>fontsize</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ad7b38f1008d32754f7f2aaa5ed8ba06e</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_RGB</type>
      <name>axis_background_col</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ac6da30ac0a3152a675cc9a8af37d91df</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_RGB &gt;</type>
      <name>data_col_list</name>
      <anchorfile>class_b_p_1_1_b_p___axes.html</anchorfile>
      <anchor>ac5555983394e879299ae996e3e45eb99</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Plot</name>
    <filename>class_b_p_1_1_b_p___plot.html</filename>
    <base>BP_Vec&lt; BP_Axes &gt;</base>
    <member kind="function">
      <type></type>
      <name>BP_Plot</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>a7b5f4753e5f67e366f370e80c492b097</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>ac190a768ce59f52a71b5031570cf0dca</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_graphic_size</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>a0d8ab57f3194d3f40509178a2752f573</anchor>
      <arglist>(double x, double y)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_graphic_margin</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>a1259e908c229dbb3f55069608774ad30</anchor>
      <arglist>(double m)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>af1c5e2815281cee291cef4ef656687c0</anchor>
      <arglist>(string filename)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_xlabel</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>a090c445d46836c14c735624df6fdf617</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_ylabel</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>a8007414b542d4007ab643f417c4cea02</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>graphic_size</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>aa417ea0eac4c14ae4612686fd0e517ee</anchor>
      <arglist>[2]</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_Text &gt;</type>
      <name>text_list</name>
      <anchorfile>class_b_p_1_1_b_p___plot.html</anchorfile>
      <anchor>a2695a5aedcd2b64a05c537d1ccfce8af</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_RVG_base</name>
    <filename>class_b_p_1_1_b_p___r_v_g__base.html</filename>
    <templarg></templarg>
    <base>BP::BP_Group</base>
    <member kind="function" virtualness="pure">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>acb39c00b44d7e74ad93348909dc99cad</anchor>
      <arglist>(BP_Gen_GVec_Member *, double)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>add_with_update</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>ae58624cfd9b08848d7dc13285388549c</anchor>
      <arglist>(BP_Gen_GVec_Member *, double)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a1ec8badea42a16572619b106a020d534</anchor>
      <arglist>(unsigned long int)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a82c3c536bd7a5ba4f6ec3b1e97b68219</anchor>
      <arglist>(BP_Gen_GVec_Member *)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>pick</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a39790848efc26cd661dd3ee8edf48c8f</anchor>
      <arglist>(double i1)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a17082811ea8ed057c7a0687d1bafaa51</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>min_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a3115f2390e3ffd1dbe49b0eff445d02c</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>max_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a686706c48a3991d8da0b3287aca0102d</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>min_rate_member</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>aff2f420f58bc687908f946e52598669f</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual BP_GVec_Member&lt; T &gt; *</type>
      <name>max_rate_member</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>aba47ed841b2cbe71af5fdf41b4c4bae1</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>refresh_total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a9fdaba878067d8d563d4436d717bf50a</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a7318f2c2a62cad73b256b432e2142df1</anchor>
      <arglist>() const =0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a4705712b33b16098f99a4e7e0489792f</anchor>
      <arglist>(unsigned long int)=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>print</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a36c69129e89024138350f0d85249f5f8</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>ad46734a97d22262e58e497824aed4b10</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a8505b778e2c87c7bea86a2759c4121f6</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a959d3494e66b4a025c4d210ce04d7e6c</anchor>
      <arglist>()=0</arglist>
    </member>
    <member kind="function" virtualness="pure">
      <type>virtual double</type>
      <name>get_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>aa32860fa960804c72540126ec760775e</anchor>
      <arglist>(unsigned long int)=0</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_RVG_base</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>a82e05ee1160b49a4faed2d0307432b87</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" virtualness="virtual">
      <type>virtual</type>
      <name>~BP_RVG_base</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__base.html</anchorfile>
      <anchor>aba74b3c2284b1087e3b7544c05d2192b</anchor>
      <arglist>()</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_RVG_linear</name>
    <filename>class_b_p_1_1_b_p___r_v_g__linear.html</filename>
    <templarg></templarg>
    <base>BP::BP_RVG_base</base>
    <member kind="function">
      <type></type>
      <name>BP_RVG_linear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a0ba45745ceffc39aa5edbac4cbc332ef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a9624b7102a1e7da5693013bc985102fc</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a35aa9caf095f83ed0c516c82c7bf2f82</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add_with_update</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a28835cbce8e53e0c3d850093e42f642c</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>ae4e6c20a952d34ea14efe976cbc83303</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a5d31fe0a0d58e286ff285f3697f6088a</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a41152980660d47db60c7cc204a9d1ebe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>ad7c9fa4c3196e4916d807f390754c680</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a3c75ffc24e3f5f2b0189267bbd53e151</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>pick</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>aa8c2669827f2e80107f7a3c0d424fc6f</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a2333506adc81bb9ee246c29a7cc25131</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>curr_total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a700ddfe0fdfdc86280c29ae153c0ed1c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a5f7c18f6f6fc7d985ec1d75fde7d5671</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>ad3892cd521882b39f233443c88b8e8b2</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>min_rate_member</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a4ce580a312ffbf2546f3bf8a1149d5c1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>max_rate_member</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a51b0f436e4c54e756da20eacd986dfae</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refresh_total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a1f901326d251e4de8d2d6707b3c29511</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>acba2d9d7ae17c613abc345c2d5cf923b</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a34873a38505b2f0bbb82960b010f77bd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a1671dcee18a9e2f93c7ded823db6c74a</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; double &gt;</type>
      <name>rates</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>ae9980a14793b98a2ea55018a10623988</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>rates_total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a2cd79ec71022b8a653a04ffcdf876de1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>rates_total_ready</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__linear.html</anchorfile>
      <anchor>a8fdf3fe2d8c0ced83d64a3247679cd91</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_RVG_tree</name>
    <filename>class_b_p_1_1_b_p___r_v_g__tree.html</filename>
    <templarg></templarg>
    <base>BP::BP_RVG_base</base>
    <member kind="function">
      <type></type>
      <name>BP_RVG_tree</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>ad8f23395b56a26fda8fe2f80e8024dc3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a8740a7551a1433fe5eade454cc283c54</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>print</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>acf13069dac41a7bbcef525df6bffaa04</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a07a393fcb8b97eca33863f6fca117d8c</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>acde7bfa3d8358d6c2dee1e2f4fa4e537</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>add_with_update</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a91f17817680175640a09c3a1194290d8</anchor>
      <arglist>(BP_Gen_GVec_Member *i1, double i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a9619395d0b651542b3124d7ffa8a90be</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a30f06dc0136c1c4239a3f7a7e296bd6d</anchor>
      <arglist>(BP_Gen_GVec_Member *i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a48ee6abd5b0b6561374681be5dd5935d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>crazy_clear</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a26263c9a3053ccd6f6a22950d6e84932</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a1bad106943d09bafaf1c25444dc9f3dd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>pick</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a2e84c5c683ee52961026cd75e4c6d8c6</anchor>
      <arglist>(double i1)</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>aaa7adf2d7e95ad63d5419a7bd9433e88</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>curr_total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a157d90003d1a555cc45f5db921552339</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>min_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a2d40ca78162044d8e59d3d1e5df03c3c</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>max_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a107cee7ccd77108439afc26fa6a306f9</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>min_rate_member</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a503eb7be761037344cf77a64a56f0846</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_GVec_Member&lt; T &gt; *</type>
      <name>max_rate_member</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a5e48f657d47c7016d7406b2eb3f200a5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>refresh_total</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>afb4aeb9ce4d5f34ed5b593441770c693</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a665b0a2c4e23633538a2dd2c3ecb54b8</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>increase_max_mvs</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a7bcc57a9ab6f36306861ecf0fbd0362e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>add_new_rate</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a7b8f4948e74fe2630888806292c76306</anchor>
      <arglist>(unsigned long int j, double r)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>clear_to_zero</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a4a7895460f166eebeeecdc5a7aab9f30</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_Vec&lt; double &gt; &gt;</type>
      <name>rate_sums</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a09812ab427a27dacc02422d8d2f49db2</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>rate_sums_ready</name>
      <anchorfile>class_b_p_1_1_b_p___r_v_g__tree.html</anchorfile>
      <anchor>a5aeaf5063bc778d71a38a2eb7403b181</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_StopWatch</name>
    <filename>class_b_p_1_1_b_p___stop_watch.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_StopWatch</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a21b10d5c59d4e3cab3eb26e760f3286c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>ab3c1689d5dae0111feb5df67b7478744</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_start</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a3b8e64ea7cf6e01cfd043335c8947865</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lap</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a2b12c99e4bf3d7c9639fc3bee480b313</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>gettime</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>abd121397ca67e4bea0f46801ea8bdec0</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>gettime_ms</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>aec61fb30013b586767b20bebd0c4bd97</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>gettime_us</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>abcf43389c1838f4c54bad1785a59c4d1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_us</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a19127f0b54a17672bffb6a6965ec089d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>get_s</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>aec608133de3ba0481489640ac987b108</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>total_time_s</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a5ab7248bc057ac0eea55573290b0d11b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>lap_time_s</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>ace36a1cec8c7211eacf48ccd201e5ad2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>start</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a5d7832af1209e6c7d970afc5d5bc1e53</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>lap</name>
      <anchorfile>class_b_p_1_1_b_p___stop_watch.html</anchorfile>
      <anchor>a67793df989ef4c6c9bb024110defd34d</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_POSCAR_class</name>
    <filename>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</filename>
    <member kind="function">
      <type>void</type>
      <name>read_POSCAR</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a3fde2cdb94f2c05bec56ad7dfa836c8c</anchor>
      <arglist>(string poscar)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_POSCAR</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac0459f3523d5a4c94d7f460f5f5c80b0</anchor>
      <arglist>(ostream &amp;sout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_POSCAR</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a4d792c1a7117b8605ebb313c77f7e04f</anchor>
      <arglist>(string poscar)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_POSCAR_class</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a1c12c81370f3482b65c377f5c99e477c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_POSCAR_class</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ae1b7f97e19e890e030e6bd85fd936977</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a15b53a94e38e27742359e66248401611</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_format_version</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a483829c3fd422c8aaf85212b7d10cf21</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_write_all_atom_types</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a55cc846cc80f251150adf34b3ae3f71e</anchor>
      <arglist>(bool b)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_lat_vectors</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ae798350b67f3c92088e3760824bb242b</anchor>
      <arglist>(double i1, const BP_Vec&lt; cart_coord_class &gt; &amp;i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>get_lat_vectors</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>aa795608f376bba3f020a546a8b402218</anchor>
      <arglist>(double &amp;d, BP_Vec&lt; cart_coord_class &gt; &amp;vec)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; cart_coord_class &gt;</type>
      <name>get_lat_vectors</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ab3d89f83b2983a3f88839a8a0fa5a05c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>get_header_line</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ae19a5bdc2078e9a192bf1456c30cbb13</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_header_line</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ae7bf51e59b53388fb1932107a0b0350f</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a97fc58ceff5819db6eba2d48c1409fb7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a528a3fe41ce370664db74270325ecbef</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a530d317d0c3de159e121c81918108b8d</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a65e1f788d383b6529d52c6f6c67eff48</anchor>
      <arglist>(string type)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>get_atom_type</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>aabaf18b8e78734155aea97d9b435cca5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>get_atom_type</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>afdb06539734e1df4e3a879444089eccf</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>atom_type_size</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>abca2c939db2a056ac51f529e2cc9710b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_atom_types_exist</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>afb972ef5115cc0a76fd5df2538acd394</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_atom_type</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ab1b5516f51cce18b185e7ff088fd0941</anchor>
      <arglist>(BP_Vec&lt; string &gt; s_list)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_atom_type</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a347b2a4e7f08a1574d6d61c8b4cadc9d</anchor>
      <arglist>(unsigned long int i, string s)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>get_type_index</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>aea43f63f7bed953083a2c1a3c063f8c8</anchor>
      <arglist>(unsigned long int i)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>get_type</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>af416ec8ab90abfd5ac47b95da565abb6</anchor>
      <arglist>(unsigned long int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_selective_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ae64a13bbda7f4595b2d49151dbd9381f</anchor>
      <arglist>(bool i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_sel_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a8374f4ed36dba91c679ac9857a82a53d</anchor>
      <arglist>(int i, BP_Vec&lt; bool &gt; sel)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_sel_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>acff4c9bf020cc12efb612c4939532521</anchor>
      <arglist>(BP_Vec&lt; BP_Vec&lt; bool &gt; &gt; sel)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_selective_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a7e85a0ff9d10fb44eedde38a3a69d040</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; BP_Vec&lt; bool &gt; &gt;</type>
      <name>get_sel_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a033a33836f746b1e70b0786a9c1a9cbd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; bool &gt;</type>
      <name>get_sel_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a3100ece866b307f2852342e01f7f7a1a</anchor>
      <arglist>(unsigned long int i)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>get_coord_mode</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a164cd87afcc72be3b9d86bf318b1115f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_coord_mode</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac3317292b4b20f44dc782475bbacd152</anchor>
      <arglist>(bool i1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_direct</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac15ce4fd49710d1cdd6fd9c55aa8f453</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>set_cart</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac85ece743933e9843a9bc94189d72c5d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_direct</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a29d05b241051ac94f0ffc4fd66564edf</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>is_cart</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ad796c8170ffa8852554fa794a325e317</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>accb9003b60b4be0305227d0474edf9e0</anchor>
      <arglist>(string type, frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a8bd14233653077dd4d95fb8e51fa68ed</anchor>
      <arglist>(unsigned long int type_index, frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac01dbe717a3066f8a5dbed295d06cc03</anchor>
      <arglist>(string type, frac_coord_class f, BP_Vec&lt; bool &gt; sel)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a70ee4e443e3c072d19b715f8da7b57fb</anchor>
      <arglist>(unsigned long int type_index, frac_coord_class f, BP_Vec&lt; bool &gt; sel)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a98c1d98da9ecde9615ee48e65c083288</anchor>
      <arglist>(string type, cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a93c995b429e3f5f3306d3f51563edc62</anchor>
      <arglist>(unsigned long int type_index, cart_coord_class c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a4af032e557cfc18e83eb4e96070077b2</anchor>
      <arglist>(string type, cart_coord_class c, BP_Vec&lt; bool &gt; sel)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a294f3818664757e8de4e03bf973a485c</anchor>
      <arglist>(unsigned long int type_index, cart_coord_class c, BP_Vec&lt; bool &gt; sel)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a5a5336fc6681e50d9c4db9310d4c3769</anchor>
      <arglist>(unsigned long int index)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove_all</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a0a0527c4a85d647b373224c1aed983f7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>in_cell</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a55ce631d8956b33f612b0190185716c9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>translate_frac</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac1adc4e2917f5862397f45b04756e907</anchor>
      <arglist>(frac_coord_class f)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>translate_cart</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a5659b8d78bbb49ae0d313636cba383ae</anchor>
      <arglist>(cart_coord_class c)</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>clean_atom_types</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ab8b780de55d66e6037f3a8bf338cc547</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function" protection="private">
      <type>void</type>
      <name>set_cum_num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ad467940a15798f91154d4e5df098c89e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>double</type>
      <name>scl</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a6754f26f4c7c576c24e1e82f99b198de</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; cart_coord_class &gt;</type>
      <name>lat_vectors</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>af7674e093cfb3a7f86e3768bdc4aaf65</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; frac_coord_class &gt;</type>
      <name>atom_pos_frac</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ad82ff614c2a9071a1b3a4e9c54df6c4f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; cart_coord_class &gt;</type>
      <name>atom_pos_cart</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ab56b60ae4ef09bab67cb3b7e495ea620</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>sum_num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ad1312bc10dcfb034f762b002e00860f8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>atom_types_exist</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>aa68bcede152f652e1b3ca670b4962ce3</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; string &gt;</type>
      <name>atom_type</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ac8a8be2e19e7f548c09c519ff30421b1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; int &gt;</type>
      <name>num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ad9aac0bf786ea57b57fb3e8eacab6101</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; int &gt;</type>
      <name>cum_num_atoms</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a27f7f59eb25a121023dde757d09f7a60</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; BP_Vec&lt; bool &gt; &gt;</type>
      <name>sel_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a8d600051f8a9fcbeab3e6ba583880736</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>direct</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>afd32ad1b494b8b08e5eb2fd6167774e8</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>header_line</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a84b75a9047d4ba489c17999469b6e611</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>bool</type>
      <name>selective_dynamics</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>a27328a67d4e0b4fd65680b5a9aac624a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>format_version</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>ae030b3eab89b81eaa5e122799bf3ddb4</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>int</type>
      <name>write_all_atom_types</name>
      <anchorfile>class_b_p_1_1_b_p___p_o_s_c_a_r__class.html</anchorfile>
      <anchor>abc77d38db1a32cce19881184401e28a7</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_Vec</name>
    <filename>class_b_p_1_1_b_p___vec.html</filename>
    <templarg>T</templarg>
    <member kind="function">
      <type>BP_Vec &amp;</type>
      <name>operator=</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a84d2e65a7a3b3793c97f603a40d56518</anchor>
      <arglist>(const BP_Vec &amp;)</arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a69b00a550cb58bcbacd79abc16ee897a</anchor>
      <arglist>(unsigned long int i1)</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>operator[]</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a19ef819840cef17fc6d59a5da53c97f1</anchor>
      <arglist>(unsigned long int i1) const </arglist>
    </member>
    <member kind="function">
      <type>T &amp;</type>
      <name>last</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>aa639907c7e88231f17031bc9fd069f95</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>const T &amp;</type>
      <name>last</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a2505527b534304efe1fbd6a4d762cf76</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>capacity</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a6d7ddcc7a1f2b70e3e5b1eb991d72e9d</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>min_capacity</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a82defa565727d822b7a935740343521f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>cap_increment</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a7be44fd5d60ac679f0639cecfeec3a05</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a41acc9f9b62cf839350c56eca4772fcd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>ab4c336a41611c5d690c9627b108f2252</anchor>
      <arglist>(const T &amp;)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>add</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>aa2ad55ebcb614b135a7a7417574ad3cf</anchor>
      <arglist>(const U &amp;)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>add_in_place</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>adb5fe5eb39621739027f1c484baf1e70</anchor>
      <arglist>(const U &amp;, unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>ad4f4cff1193bef9ffc018c4b5c8fb310</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ordered_remove</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>ad2e87024a9b6e34537ed8aeea6beed37</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>swap</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a23107e3984c2a61706b5d61be53b7270</anchor>
      <arglist>(unsigned long int, unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>ordered_swap</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a65e20407be59c1ac70b85b8d4209cea8</anchor>
      <arglist>(unsigned long int i1, unsigned long int i2)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>ada83d75d4f3ae5610b3d41d7d10033a3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>free</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a938d34802776d288f2af45176295a63e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>erase</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a7a1f33553fd300bbebb2946e6c094efd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a87835d8a2c685bdb758102d94fb14e54</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size_alloc</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>ac81ced9e8899c2ea2941b3ac19dfcec0</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>size_capacity</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>ad2a4f090804688402a40ae1ace2dde22</anchor>
      <arglist>() const </arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>put_in</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>afd4157fd49aaf9dd5ef2763e4f9b7e30</anchor>
      <arglist>(T *)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>put_in_place</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a79b039b5ebc186ed99c0c7026a0ab7f5</anchor>
      <arglist>(T *, unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>take_out</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>afb516b23909dfba287a51ca05962cf5c</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>T *</type>
      <name>ordered_take_out</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a368bda1490ad7b2b9e327f49fe24c370</anchor>
      <arglist>(unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec</type>
      <name>range</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>af791d60d533dcefb60375dcd93669422</anchor>
      <arglist>(unsigned long int, unsigned long int)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>append</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>aed2bb30d526e1c48b8d06bc352989b5d</anchor>
      <arglist>(const BP_Vec &amp;)</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>find_first</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a9d89cc8e4f80785f5423e1705990caab</anchor>
      <arglist>(const T &amp;)</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>find_first</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a572602bbace09c71e0f8636f4a1599f9</anchor>
      <arglist>(const T &amp;i1, unsigned long int &amp;index)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; unsigned long int &gt;</type>
      <name>find_all</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>aea60f9eaba064f4780692777e848dfd2</anchor>
      <arglist>(const T &amp;)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_multiline</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>aa49acf439e9a0723b4e6656b40755978</anchor>
      <arglist>(ostream &amp;sout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_inline</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a6c8d995fd5368165519f70f62c424322</anchor>
      <arglist>(ostream &amp;sout)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>init</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a1a2d82ebf62aab28c4cdef8dc0bf116b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Vec</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a0fd97358dc5cd538b7cfc9b8d553c893</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Vec</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a050671240a283d02c7014bab38ab7db9</anchor>
      <arglist>(const BP_Vec &amp;)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_Vec</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a6ffbfd03dcb837eba817fde2334b8d9f</anchor>
      <arglist>(unsigned long int i1, const T &amp;i2)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_Vec</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>af6ec1e265aeed529bfa792516544f587</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>N_max</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>adcf81e910aa995350f44f2a2c5f3aaf1</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>N_alloc</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a5438e97e98197fcdbe0e81ccdaa61a6a</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>unsigned long int</type>
      <name>N</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a72b40dcea72d2df2aabd930d16d4976f</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>int</type>
      <name>cap_incr</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a98b617decee6c9a62972d1601b62e49d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="protected">
      <type>T **</type>
      <name>val</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>af2a812809f195ba3fce75d04b7c7daf8</anchor>
      <arglist></arglist>
    </member>
    <member kind="friend">
      <type>friend ostream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1_b_p___vec.html</anchorfile>
      <anchor>a5ea4534e952534f5613b973494755f6e</anchor>
      <arglist>(ostream &amp;outstream, const BP_Vec&lt; U &gt; &amp;vec)</arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_zParse</name>
    <filename>class_b_p_1_1_b_p__z_parse.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_zParse</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>ab24936d5a06f0e085e7749321fde1998</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_zParse</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a401b0cd1e9c967946a7965e5585ae559</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>aa8f41d744e0317cd811a7c64c8d5138c</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a391f18343ea4157d3056a9d055edf7e8</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>add_com</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>ac15d8d3646b03885f929435100c405b0</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>clear_com</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a7368a4ab34d3257d2b86e37a8d29218e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>remove_com</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a712049f3bd7aaae6681ba23663741575</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>get_com</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a4df1e7ef969f4e4bce5c4264039268d7</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>com</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a08423ab1a8fbf14ef2dc25194bd6c16a</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getline</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a546e69605a7cb827fcb5661e079737dd</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>getline_all</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>aaf304143429bf396e0d81d47c3ffef36</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; string &gt;</type>
      <name>getline_string</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a8495ff0fead9e2f63622ae39f30a04a9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; int &gt;</type>
      <name>getline_int</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a5d752371dac11f2f37afe9bbb9c19696</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>BP_Vec&lt; double &gt;</type>
      <name>getline_double</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a15ad136933d1cf885b5dfc64818deac9</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>next_string</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a16c8dfa2022f852cf085de0bcd68f1d3</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>next_int</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a685382b39af1fe0f1fdf606a78ccb28b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>long int</type>
      <name>next_lint</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>abc3426d714042f5da324ce2b0cbe3388</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>next_ulint</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a17b720ad697eafc1da3775110b509ad2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>next_double</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>adc71f0862ba18785cd65e48c92fd3538</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>eof</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a2b86d05ab13c47181dfd1f7a7ee4c69e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>igzstream &amp;</type>
      <name>operator&gt;&gt;</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>ac413a366284b128319e116c3559c3846</anchor>
      <arglist>(T &amp;t)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>igzstream</type>
      <name>infile</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>ac0a8cac88b89b2d9ed13417d3b34795d</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>BP_Vec&lt; string &gt;</type>
      <name>com_list</name>
      <anchorfile>class_b_p_1_1_b_p__z_parse.html</anchorfile>
      <anchor>a1ce025d807ed99459b6f5615afe36d10</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_bzParse</name>
    <filename>class_b_p_1_1_b_p__bz_parse.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_bzParse</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>abc4fd028ed10a86593fc9ab6b2ab41a6</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_bzParse</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>afdac7c7743d418056a2b053b96eec98a</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>a56162f31e1a4520d3f5ee67938f2de9d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>reset</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>ad3d2b1458d8505b71d42ae5125ea07a7</anchor>
      <arglist>(string s1)</arglist>
    </member>
    <member kind="function">
      <type>char</type>
      <name>next_char</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>ae60f06c15bd603138600e835876aae65</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>int</type>
      <name>next_int</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>a5d15b01a0a9d7135cf2b1ecc6d909531</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>unsigned long int</type>
      <name>next_ulint</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>adfaf6a175b9efa99e3f964ebbb9ef9e5</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>long int</type>
      <name>next_lint</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>a009e3aa99c29d16abb7bdfda8ea3a999</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>double</type>
      <name>next_double</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>a9d64f89465bdb79a7b9ee2d0ca18d5e1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>peek</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>a8237abe8f6ab78af14d1f0d0d2c5c382</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>bool</type>
      <name>eof</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>ab19c2165d1aea0ddeb3410f26898cf3f</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>igzstream</type>
      <name>infile</name>
      <anchorfile>class_b_p_1_1_b_p__bz_parse.html</anchorfile>
      <anchor>ac88d4b79bfa452a2f389b4dec7877032</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_zWrite</name>
    <filename>class_b_p_1_1_b_p__z_write.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_zWrite</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>af922c1d21d9faf39a9542f01f30fa31b</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_zWrite</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a99477e6430543c89f320e98eeab82a3c</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_zWrite</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a417822382e217513dc759139ee08a365</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a90243403e6642044994090c89a415748</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a6b36b390dd36e71f1ac683486b3227fa</anchor>
      <arglist>(string s, string mode)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>newfile</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a57725b9adcef72598e29c707a17ed5e1</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>name</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a7ef74ad4f3dec4dea647200a009aa4fe</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>afbc734d4951771241855968b32f2acce</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>ogzstream &amp;</type>
      <name>operator&lt;&lt;</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>abd57413e2179c2b65a3771ac0d9bb7bc</anchor>
      <arglist>(const T &amp;t)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ogzstream</type>
      <name>outfile</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>aa4d0dcd2ab5506513503c2ced0bbc224</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>filename</name>
      <anchorfile>class_b_p_1_1_b_p__z_write.html</anchorfile>
      <anchor>a5c0f20af7d50978381c40cf3e3367bbd</anchor>
      <arglist></arglist>
    </member>
  </compound>
  <compound kind="class">
    <name>BP::BP_bzWrite</name>
    <filename>class_b_p_1_1_b_p__bz_write.html</filename>
    <member kind="function">
      <type></type>
      <name>BP_bzWrite</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>abd1bba006323077d5d8e30cf739c67e2</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>BP_bzWrite</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>ad3d9d8888e5c205d31fe1548e251f5ce</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type></type>
      <name>~BP_bzWrite</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>ab94faad0596abe0d482c298ad097fb9e</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a92eaadea5f0f8ee07a34d404291f2d2a</anchor>
      <arglist>(string s)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>changefile</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>af9c15d0079a21c354993f286bbf1e39d</anchor>
      <arglist>(string s, string mode)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>newfile</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a0834079adb9818189f6dfafb257c4c45</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>string</type>
      <name>name</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>acfe01bce4c428fed4688147865c15e00</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>close</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a6af0f6e27544d64e4e06d935aa811c8d</anchor>
      <arglist>()</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>ae4baee679b2975dae73537f327f6d286</anchor>
      <arglist>(const char *s, unsigned long int n)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_char</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a7e997b65d93de2ff5740a1bfcffab754</anchor>
      <arglist>(char c)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_int</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a098af22743d6a953c4c8442d7e0918d3</anchor>
      <arglist>(int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_lint</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a158f0d2531faaa70f81d557a78be7a7e</anchor>
      <arglist>(long int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_ulint</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a2b99c086d54554199db8510db286748b</anchor>
      <arglist>(unsigned long int i)</arglist>
    </member>
    <member kind="function">
      <type>void</type>
      <name>write_double</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a68ee5d17b2e89daf5d9707cfa9ed0e81</anchor>
      <arglist>(double d)</arglist>
    </member>
    <member kind="variable" protection="private">
      <type>ogzstream</type>
      <name>outfile</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a9a24c945bbb0e6068b4c8c056a1a945c</anchor>
      <arglist></arglist>
    </member>
    <member kind="variable" protection="private">
      <type>string</type>
      <name>filename</name>
      <anchorfile>class_b_p_1_1_b_p__bz_write.html</anchorfile>
      <anchor>a469054393bc986ab19030bb7fa523c72</anchor>
      <arglist></arglist>
    </member>
  </compound>
</tagfile>
