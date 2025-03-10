import { defineConfig } from 'vitepress'
import { tabsMarkdownPlugin } from 'vitepress-plugin-tabs'
import mathjax3 from "markdown-it-mathjax3";
import footnote from "markdown-it-footnote";
import path from 'path'

function getBaseRepository(base: string): string {
  if (!base || base === '/') return '/';
  const parts = base.split('/').filter(Boolean);
  return parts.length > 0 ? `/${parts[0]}/` : '/';
}

const baseTemp = {
  base: '/MixedModels.jl/',// TODO: replace this in makedocs!
}

const navTemp = {
  nav: [
{ text: 'MixedModels.jl Documentation', link: '/index' },
{ text: 'Articles', collapsed: false, items: [
{ text: 'Model constructors', link: '/constructors' },
{ text: 'Details of the parameter estimation', link: '/optimization' },
{ text: 'Normalized Gauss-Hermite Quadrature', link: '/GaussHermite' },
{ text: 'Prediction and simulation in Mixed-Effects Models', link: '/prediction' },
{ text: 'Parametric bootstrap for mixed-effects models', link: '/bootstrap' },
{ text: 'Rank deficiency in mixed-effects models', link: '/rankdeficiency' },
{ text: 'Alternative display and output formats', link: '/mime' }]
 },
{ text: 'API', link: '/api' }
]
,
}

const nav = [
  ...navTemp.nav,
  {
    component: 'VersionPicker'
  }
]

// https://vitepress.dev/reference/site-config
export default defineConfig({
  base: '/MixedModels.jl/',// TODO: replace this in makedocs!
  title: 'MixedModels',
  description: 'Documentation for MixedModels.jl',
  lastUpdated: true,
  cleanUrls: true,
  outDir: '../final_site', // This is required for MarkdownVitepress to work correctly...
  head: [
    
    ['script', {src: `${getBaseRepository(baseTemp.base)}versions.js`}],
    // ['script', {src: '/versions.js'], for custom domains, I guess if deploy_url is available.
    ['script', {src: `${baseTemp.base}siteinfo.js`}]
  ],
  
  vite: {
    resolve: {
      alias: {
        '@': path.resolve(__dirname, '../components')
      }
    },
    optimizeDeps: {
      exclude: [ 
        '@nolebase/vitepress-plugin-enhanced-readabilities/client',
        'vitepress',
        '@nolebase/ui',
      ], 
    }, 
    ssr: { 
      noExternal: [ 
        // If there are other packages that need to be processed by Vite, you can add them here.
        '@nolebase/vitepress-plugin-enhanced-readabilities',
        '@nolebase/ui',
      ], 
    },
  },
  markdown: {
    math: true,
    config(md) {
      md.use(tabsMarkdownPlugin),
      md.use(mathjax3),
      md.use(footnote)
    },
    theme: {
      light: "github-light",
      dark: "github-dark"}
  },
  themeConfig: {
    outline: 'deep',
    logo: { src: '/logo.png', width: 24, height: 24},
    search: {
      provider: 'local',
      options: {
        detailedView: true
      }
    },
    nav,
    sidebar: [
{ text: 'MixedModels.jl Documentation', link: '/index' },
{ text: 'Articles', collapsed: false, items: [
{ text: 'Model constructors', link: '/constructors' },
{ text: 'Details of the parameter estimation', link: '/optimization' },
{ text: 'Normalized Gauss-Hermite Quadrature', link: '/GaussHermite' },
{ text: 'Prediction and simulation in Mixed-Effects Models', link: '/prediction' },
{ text: 'Parametric bootstrap for mixed-effects models', link: '/bootstrap' },
{ text: 'Rank deficiency in mixed-effects models', link: '/rankdeficiency' },
{ text: 'Alternative display and output formats', link: '/mime' }]
 },
{ text: 'API', link: '/api' }
]
,
    editLink: { pattern: "https://https://github.com/ajinkya-k/MixedModels.jl/edit/main/docs/src/:path" },
    socialLinks: [
      { icon: 'github', link: 'https://github.com/ajinkya-k/MixedModels.jl' }
    ],
    footer: {
      message: 'Made with <a href="https://luxdl.github.io/DocumenterVitepress.jl/dev/" target="_blank"><strong>DocumenterVitepress.jl</strong></a><br>',
      copyright: `© Copyright ${new Date().getUTCFullYear()}.`
    }
  }
})
