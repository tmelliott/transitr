/* Mock Realtime data server
 *
 * Serve "realtime" GTFS data by request? i.e., use worker's pid to serve the next file for that worker.
 */

const express = require('express')
const app = express()

// load all of the files in the ZIP archieve ...
const fs = require('fs')
const path = require('path')
var all_files = fs.readdirSync('archive')
const files = all_files.filter(file => file.match('^vehicle'))
const files2 = all_files.filter(file => file.match('^trip'))

var pids = [] // {"pid001": 1}

app.use('/archive', express.static('archive'))

app.use("/:pid/:ts/minutes/:min/vehicle_positions", (req, res) => {
    const pid = req.params.pid
    var ts = parseInt(req.params.ts)
    const min = parseInt(req.params.min)
    ts -= min*60
    ts *= 1000
    var x = new Date()
    x.setTime(ts)
    const d = x.getFullYear() + 
        (x.getMonth() < 9 ? 0 : '') + (x.getMonth() + 1) + 
        (x.getDate() < 10 ? 0 : '') + x.getDate() +
        (x.getHours() < 10 ? 0 : '' ) + x.getHours() +
        (x.getMinutes() < 10 ? 0 : '' ) + x.getMinutes() +
        (x.getSeconds() < 10 ? 0 : '' ) + x.getSeconds()
    const dt = parseInt(d)
    var start = 0
    for (var i=0; i<files.length; i++) {
        if (parseInt(files[i].replace(/[a-z_\.]*/g, '')) > dt) {
            start = i-1
            break
        }
    }
    if (!(pid in pids)) pids[pid] = 0
    if (pids[pid] < start) {
        console.log('starting at ' + start)
        pids[pid] = start
    }
    if (pids[pid] >= files.length) {
        res.send('simulation complete')
        return
    }
    res.download(path.join(__dirname, 'archive/' + files[pids[pid]]));
    pids[pid]++
})

app.use('/:pid/:start/:end/vehicle_positions', (req, res) => {
    const pid = req.params.pid
    const start = parseInt(req.params.start)
    const end = parseInt(req.params.end)
    if (!(pid in pids)) pids[pid] = start
    if (pids[pid] < start) {
        pids[pid] = start
    }
    if (pids[pid] >= end) {
        res.send('simulation complete')
        return
    }
    res.download(path.join(__dirname, 'archive/' + files[pids[pid]]));
    pids[pid]++
})

app.use('/:pid/vehicle_positions', (req, res) => {
    const pid = req.params.pid
    if (!(pid in pids)) pids[pid] = 0
    if (pids[pid] >= files.length) {
        res.send('simulation complete')
        return
    }
    res.download(path.join(__dirname, 'archive/' + files[pids[pid]]));
    pids[pid]++

    // console.log(pids)
})

app.use('/:pid/trip_updates', (req, res) => {
    const pid = req.params.pid
    if (!(pid in pids)) pids[pid] = 1
    if (pids[pid] - 1 >= files.length) {
        res.send('simulation complete')
        return
    }
    res.download(path.join(__dirname, 'archive/' + files2[pids[pid]-1]));
})

app.use('/:pid/reset', (req, res) => {
    const pid = req.params.pid
    if (pid in pids) pids[pid] = 0
    res.send('ok')
})

app.get('*', (req, res) => {
    console.log(req)
    res.status(404).send('Not found')
})

app.listen(3000, () => console.log('Mock GTFS server running on port 3000!'))
