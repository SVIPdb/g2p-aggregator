# -*- coding: utf-8 -*-

from random import choice
import string

from vdom import style, div, a, h1, img

def expando(summary, expanded, content):
    prefix = ''.join([choice(string.ascii_uppercase + string.digits) for _ in range(16)])
    
    noAStyle = dict(textDecoration='none', color="currentColor")

    return div(
        a(summary, href="#%s-summary" % prefix, id="%s-summary" % prefix, className="hidden", style=noAStyle),
        a(expanded, href="#%s-c" % prefix, id="%s-c" % prefix, className="shown", style=noAStyle),
        div(content, className="content"),
        style("""
        .content, .shown, .hidden:target {
          display: none;
        }
        .hidden:target + .shown, .hidden:target ~ .content {
          display: block;
          margin-left: 10px;
        }
        """
        )
    )

def collapse_xml(tree):
    if len(tree):
        return div(**(expando(u"▶ %s" % elem.name, u"▼ %s" % elem.name, collapse_xml(elem)) for elem in tree))
    else:
        return div(tree.name)
