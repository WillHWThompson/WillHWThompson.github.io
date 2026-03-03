---
layout: post-fullwidth
title: "Scientific Workflows on the Eve of the Butlerian Jihad"
date: 2026-02-28 09:00:00
description: >
  A practical manifesto for scientific software in the era of vibecoding:
  let's write beautiful code while we still can.
tags: scientific-computing reproducibility workflow python
categories: research
---

<iframe
  id="post-iframe"
  src="{{ '/assets/posts/scientific-software/index.html' | relative_url }}"
  style="width: 100%; height: 8000px; border: none; display: block; overflow: hidden;"
  scrolling="no"
  title="Scientific Workflows on the Eve of the Butlerian Jihad"
></iframe>
<script>
(function() {
  var iframe = document.getElementById('post-iframe');
  function resize() {
    try {
      var h = iframe.contentDocument.documentElement.scrollHeight;
      if (h > 100) iframe.style.height = h + 'px';
    } catch(e) {}
  }
  iframe.addEventListener('load', function() {
    resize();
    setTimeout(resize, 500);
    setTimeout(resize, 1500);
    setTimeout(resize, 4000);
  });
})();
</script>
