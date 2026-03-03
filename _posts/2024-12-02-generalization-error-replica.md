---
layout: post-fullwidth
title: "Generalization Error via the Replica Method"
date: 2024-12-02 12:00:00
description: >
  A statistical-physics treatment of high-dimensional learning using the replica
  method from disordered systems.
tags: machine-learning statistical-physics random-matrix
categories: research
thumbnail: assets/posts/generalization-error/figures/double_descent.png
---

<iframe
  id="post-iframe"
  src="{{ '/assets/posts/generalization-error/index.html' | relative_url }}"
  style="width: 100%; height: 5000px; border: none; display: block; overflow: hidden;"
  scrolling="no"
  title="Generalization Error via the Replica Method"
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
