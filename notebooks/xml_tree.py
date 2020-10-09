from lxml import etree
from IPython.display import HTML

def make_node(elem):
    kids = ('<div class="kid_holder">%s</div>' % "".join((make_node(k) for k in elem))
        if len(elem)
        else "")
    
    attrs = '<span class="attr_list">%s</span>' % (" ".join(
        (' <span class="key">%s</span>=<span class="value">"%s"</span>' % (k,v)) for k,v in elem.attrib.items()
    ) if len(elem.attrib) > 0 else "")
    
    postamble = f"""<div class="tag post" onClick="toggle_elem_contents(this)">&lt;/{elem.tag}&gt;</div>"""

    return (f"""
    <div class="node">
        <div class="tag" onClick="toggle_elem_contents(this)">&lt;{elem.tag}{attrs}{'/' if not len(elem) and not elem.text else ""}&gt;</div>
        {kids}
        {'<span class="text">%s</span>' % elem.text if elem.text and elem.text.strip() != "" else ""}
        </span>{postamble if len(elem) or elem.text else ""}
    </div>
    """)

def render_tree(root):
    return HTML("""
    <script>
    function toggle_elem_contents(elem) {
        const container = elem.parentElement;
        container.classList.toggle('hidden');
    }
    </script>
    <div class="tree">""" + make_node(root) + """</div>
    <style>
    .tree { font-family: Menlo, sans-serif; }
    .tag { cursor: pointer; font-weight: bold; color: green; display: inline-block; }
    .kid_holder { margin-left: 10px; border-left: dotted 1px gray; padding-left: 5px; }
    .text { color: #333; display: inline-block; }
    
    .attr_list { font-weight: normal; }
    .attr_list .key {color: #7D9029; }
    .attr_list .value {color: darkred; }
        
    .hidden .tag { color: gray; }
    .hidden .tag.post::before { content: '... '; }
    .hidden .kid_holder { display: none; }
    .hidden .text { display: none; }
    </style>
    """)
